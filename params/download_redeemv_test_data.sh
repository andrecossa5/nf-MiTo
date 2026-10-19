#!/usr/bin/env bash
# Download the small ReDeeM test data bundled with the RedeemV repository (commit 5c34940):
# 250,000 read triplets of a 10x Multiome ATAC library enriched for chrM (R1, 24-nt barcode read i5, R2) and the
# Cell Ranger ARC cell list (10,125 barcodes, RNA space), ~31 MB.
# Writes <outdir>/fastq/quickstart_S1_L001_R{1,2,3}_001.fastq.gz (cellranger-atac names: R2 = barcode read, R3 = second
# genomic read), <outdir>/cells_rna.txt and <outdir>/samplesheet.csv. Run nf-MiTo with --assay ReDeeM --barcode_space RNA.
# Usage: bash download_redeemv_test_data.sh [outdir]
set -euo pipefail

OUT=$(mkdir -p "${1:-redeemv_test_data}/fastq" && cd "${1:-redeemv_test_data}" && pwd)
BASE=https://raw.githubusercontent.com/chenweng1991/redeemV/5c34940/source

curl -fsSL -o "$OUT/fastq/quickstart_S1_L001_R1_001.fastq.gz" "$BASE/Example.R1.fq.gz"
curl -fsSL -o "$OUT/fastq/quickstart_S1_L001_R2_001.fastq.gz" "$BASE/Example.i5.fq.gz"
curl -fsSL -o "$OUT/fastq/quickstart_S1_L001_R3_001.fastq.gz" "$BASE/Example.R2.fq.gz"
curl -fsSL "$BASE/barcodes.tsv.gz" | gzip -dc > "$OUT/cells_rna.txt"

printf 'sample,fastq_dir,bam,cell_barcodes\nquickstart,%s/fastq,,%s/cells_rna.txt\n' "$OUT" "$OUT" > "$OUT/samplesheet.csv"
echo "Done: $OUT ($(wc -l < "$OUT/cells_rna.txt") cell barcodes)"; cat "$OUT/samplesheet.csv"
