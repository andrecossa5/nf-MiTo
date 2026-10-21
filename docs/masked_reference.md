# Making a masked reference

`mtscATAC` and `ReDeeM` align reads to a GRCh38 genome in which the **nuclear mitochondrial insertions (NUMTs)** are masked, so that nuclear copies of mtDNA do not contribute reads to the MT BAM. `MAESTER` uses the same FASTA to build its STAR index. nf-MiTo does not download or mask a genome: you give the masked FASTA with `--masked_ref`.

The FASTA must be **uncompressed**, must contain the MT contig (named as in `--string_MT`, `chrM` by default), and the masked regions are replaced by `N`.

## Requirements

[bedtools](https://bedtools.readthedocs.io/) and, optionally, `pigz`:

```bash
mamba create -n masking -c bioconda -c conda-forge bedtools pigz
mamba activate masking
```

## Steps

We use the GENCODE release 32 genome (GRCh38.p13, `chr`-prefixed names) and the NUMT intervals of the [mitoblacklist](https://github.com/caleblareau/mitoblacklist) repository (`hg38.full.blacklist.bed`, 820 nuclear intervals).

```bash
mkdir -p reference && cd reference

# 1. The genome, uncompressed
wget -qO- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.p13.genome.fa.gz \
  | pigz -dc > genome.fa

# 2. The NUMT intervals (BED)
wget -q https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed

# 3. Mask them with N
bedtools maskfasta -fi genome.fa -bed hg38.full.blacklist.bed -fo masked.fa
```

## Check

```bash
grep '>' masked.fa | grep -w chrM       # the MT contig is there, with the name of --string_MT
```

Then:

```json
{ "masked_ref": "/path/to/reference/masked.fa" }
```

The bwa-mem2 index (and, for MAESTER, the STAR index) is built by the pipeline from this FASTA, once, and stored in `<output_folder>/resources` for reuse.

## Other genomes

Any other reference works the same way, as long as you have a BED file with its NUMT intervals: change the genome URL and the BED in steps 1 and 2, and set `--string_MT` to the name of the MT contig of your FASTA. Use the same masked FASTA for MAESTER, for the STAR index (see [Building a STAR index](star_index.md)) and for every other assay of a project.
