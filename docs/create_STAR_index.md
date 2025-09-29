# Creating a STAR Index for nf-MiTo

This tutorial provides step-by-step instructions for creating a custom STAR index for use with the nf-MiTo pipeline. This is useful when you want to:

- Use a custom reference genome not available in standard repositories
- Pre-build an index to speed up multiple pipeline runs
- Use specific genome versions or annotations
- Work in environments with limited internet connectivity

## Prerequisites

Before starting, ensure you have the following tools installed:

- **STAR** (≥2.7.0) - RNA-seq aligner
- **bedtools** (≥2.26.0) - Genome arithmetic tools  
- **wget** - File download utility
- **pigz** - Parallel gzip implementation (optional but recommended)

You can install these using mamba:

```bash
# Create a new conda environment
mamba create -n star-index star bedtools wget pigz -c bioconda -c conda-forge

# Activate the environment
conda activate star-index
```

## Overview

The STAR index creation process involves three main steps:

1. **Prepare the reference genome** - Download and mask NUMT regions
2. **Prepare the GTF annotation** - Download and filter gene annotations
3. **Build the STAR index** - Generate the searchable genome index

## Step 1: Prepare the Reference Genome

### 1.1 Set up working directory

```bash
# Create working directory
mkdir -p star_index_creation
cd star_index_creation

# Create subdirectories
mkdir -p downloads processed_files STAR_index
```

### 1.2 Download the reference genome

```bash
# Default human reference genome (GRCh38.p13)
REFERENCE_GENOME="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.p13.genome.fa.gz"

# Download and decompress the genome
echo "Downloading reference genome..."
wget -qO- ${REFERENCE_GENOME} | pigz -dc > downloads/genome.fa
echo "Genome download completed: $(wc -l < downloads/genome.fa) lines"
```

### 1.3 Download and apply NUMT masking

Nuclear mitochondrial DNA sequences (NUMTs) can interfere with mitochondrial analysis, so we mask them:

```bash
# Download NUMT blacklist regions
NUMTS_REGIONS="https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed"

echo "Downloading NUMT regions..."
wget ${NUMTS_REGIONS} -O downloads/hg38.full.blacklist.bed

# Mask NUMT regions in the genome
echo "Masking NUMT regions..."
bedtools maskfasta \
  -fi downloads/genome.fa \
  -bed downloads/hg38.full.blacklist.bed \
  -fo processed_files/masked_genome.fa

echo "NUMT masking completed: $(grep -c '^>' processed_files/masked_genome.fa) sequences"
```

## Step 2: Prepare the GTF Annotation

### 2.1 Download the GTF file

```bash
# Default GTF annotation (GENCODE v32)
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz"

echo "Downloading GTF annotation..."
wget -qO- ${GTF_URL} | gunzip > downloads/genes.gtf
echo "GTF download completed: $(wc -l < downloads/genes.gtf) lines"
```

### 2.2 Filter the GTF annotation

We filter the GTF to include only relevant gene types for better performance:

```bash
# Create GTF filtering script
cat > filter_gtf.sh << 'EOF'
#!/bin/bash

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input.gtf> <output.gtf>"
  exit 1
fi

INPUT="$1"
OUTPUT="$2"

echo "Filtering GTF: $INPUT -> $OUTPUT"

# Copy header lines
grep '^#' "$INPUT" > "$OUTPUT" || true

# Filter GTF to include relevant gene types
awk -F'\t' 'BEGIN {OFS="\t"} 
!/^#/ {
    # Look for gene_type or gene_biotype in the attributes column
    if ($9 ~ /gene_type "protein_coding"/ || 
        $9 ~ /gene_type "lncRNA"/ || 
        $9 ~ /gene_type "antisense"/ ||
        $9 ~ /gene_type "IG_/ ||
        $9 ~ /gene_type "TR_/ ||
        $9 ~ /gene_biotype "protein_coding"/ || 
        $9 ~ /gene_biotype "lncRNA"/ || 
        $9 ~ /gene_biotype "antisense"/ ||
        $9 ~ /gene_biotype "IG_/ ||
        $9 ~ /gene_biotype "TR_/) {
        print $0
    }
}' "$INPUT" >> "$OUTPUT"

echo "GTF filtering completed"
EOF

# Make script executable
chmod +x filter_gtf.sh

# Filter the GTF
echo "Filtering GTF annotation..."
./filter_gtf.sh downloads/genes.gtf processed_files/filtered_genes.gtf

echo "GTF filtering completed: $(wc -l < processed_files/filtered_genes.gtf) lines retained"
```

### 2.3 Verify GTF filtering (optional)

```bash
# Check gene types before filtering
echo "=== Gene types in original GTF ==="
grep -v '^#' downloads/genes.gtf | \
  grep -E 'gene_type|gene_biotype' | \
  sed -E 's/.*gene_(type|biotype) "([^"]+)".*/\2/' | \
  sort | uniq -c | sort -nr | head -10

echo -e "\n=== Gene types in filtered GTF ==="
grep -v '^#' processed_files/filtered_genes.gtf | \
  grep -E 'gene_type|gene_biotype' | \
  sed -E 's/.*gene_(type|biotype) "([^"]+)".*/\2/' | \
  sort | uniq -c | sort -nr
```

## Step 3: Build the STAR Index

### 3.1 Determine optimal parameters

```bash
# Check available memory and CPU cores
echo "Available system resources:"
echo "Memory: $(free -h | grep '^Mem:' | awk '{print $2}')"
echo "CPU cores: $(nproc)"

# Calculate recommended genomeSAindexNbases
# For human genome: typically 14 (default is good for most genomes)
GENOME_SIZE=$(wc -c < processed_files/masked_genome.fa)
SA_INDEX_NBASES=$(python3 -c "import math; print(min(14, int(math.log2(${GENOME_SIZE})/2 - 1)))")

echo "Genome size: ${GENOME_SIZE} bytes"
echo "Recommended genomeSAindexNbases: ${SA_INDEX_NBASES}"
```

### 3.2 Build the STAR index

```bash
# Set number of threads (adjust based on your system)
THREADS=$(nproc)

echo "Building STAR index with ${THREADS} threads..."

# Build STAR index
STAR \
  --runMode genomeGenerate \
  --genomeDir STAR_index \
  --genomeFastaFiles processed_files/masked_genome.fa \
  --sjdbGTFfile processed_files/filtered_genes.gtf \
  --runThreadN ${THREADS} \
  --genomeSAindexNbases ${SA_INDEX_NBASES}

echo "STAR index build completed!"
```

### 3.3 Verify the index

```bash
# Check index files
echo "=== STAR index files created ==="
ls -la STAR_index/

# Verify essential files are present
REQUIRED_FILES=("Genome" "SA" "SAindex" "chrLength.txt" "chrName.txt" "chrStart.txt" "genomeParameters.txt")

echo -e "\n=== Verifying required files ==="
for file in "${REQUIRED_FILES[@]}"; do
  if [[ -f "STAR_index/$file" ]]; then
    echo "✓ $file - $(ls -lh STAR_index/$file | awk '{print $5}')"
  else
    echo "✗ $file - MISSING"
  fi
done

# Display index statistics
echo -e "\n=== Index statistics ==="
if [[ -f "STAR_index/Log.out" ]]; then
  grep -E "(Number of|Genome|SA" STAR_index/Log.out
fi
```

## Step 4: Using the Index with nf-MiTo

### 4.1 Copy index to desired location

```bash
# Copy to your nf-MiTo resources directory
FINAL_INDEX_PATH="/path/to/your/nf-MiTo/resources/STAR_index"

# Create directory and copy
mkdir -p "$(dirname ${FINAL_INDEX_PATH})"
cp -r STAR_index "${FINAL_INDEX_PATH}"

echo "Index copied to: ${FINAL_INDEX_PATH}"
```

### 4.2 Configure nf-MiTo parameters

Update your nf-MiTo parameter file to use the prebuilt index:

```json
{
  "build_STAR_index": false,
  "prebuilt_STAR_index": "/path/to/your/nf-MiTo/resources/STAR_index",
  "string_MT": "chrM"
}
```

Or specify via command line:

```bash
nextflow run main.nf \
  --build_STAR_index false \
  --prebuilt_STAR_index /path/to/your/nf-MiTo/resources/STAR_index \
  --raw_data_input your_samples.csv \
  --output_folder results
```

## Advanced Options

### Custom Reference Genomes

To use a different reference genome:

```bash
# Example: Using mouse genome
REFERENCE_GENOME="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.p6.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz"

# Note: You'll need a mouse-specific NUMT blacklist or skip masking
```

### Optimizing for Different Systems

**For high-memory systems (>64GB RAM):**

```bash
STAR \
  --runMode genomeGenerate \
  --genomeDir STAR_index \
  --genomeFastaFiles processed_files/masked_genome.fa \
  --sjdbGTFfile processed_files/filtered_genes.gtf \
  --runThreadN ${THREADS} \
  --limitGenomeGenerateRAM 50000000000
```

**For limited-memory systems (<16GB RAM):**

```bash
STAR \
  --runMode genomeGenerate \
  --genomeDir STAR_index \
  --genomeFastaFiles processed_files/masked_genome.fa \
  --sjdbGTFfile processed_files/filtered_genes.gtf \
  --runThreadN ${THREADS} \
  --genomeSAindexNbases 12 \
  --genomeChrBinNbits 15
```

## Summary

You now have a complete STAR index ready for use with nf-MiTo! The key files are:

- **STAR_index/** - The complete index directory
- **processed_files/masked_genome.fa** - NUMT-masked reference genome
- **processed_files/filtered_genes.gtf** - Filtered gene annotations

This index can be reused for multiple nf-MiTo runs and shared across projects using the same reference genome.

For more information about nf-MiTo parameters and workflows, see the main [README](../README.md) and [PREPROCESS tutorial](PREPROCESS.md).
