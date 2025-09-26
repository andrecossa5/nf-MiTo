#!/bin/bash

# Usage: bash filter_gtf.sh input.gtf output.gtf

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

# Filter GTF using basic AWK compatible with older versions
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