#!/usr/bin/env bash
# Download the small mtscATAC test BAM bundled with the mgatk repository (commit 7ca7d74, mgatk 0.9.0):
# cra_v2.bam, 50,000 MT reads of 37 cells (10x scATAC, GRCh38, chrM, barcode in CB), 2.7 MB.
# Writes <outdir>/{cra_v2.bam,cra_v2.bam.bai,cra_v2_barcodes.txt,samplesheet.csv}.
# Usage: bash download_mgatk_test_data.sh [outdir]
set -euo pipefail

OUT=$(mkdir -p "${1:-mgatk_test_data}" && cd "${1:-mgatk_test_data}" && pwd)
BASE=https://raw.githubusercontent.com/caleblareau/mgatk/7ca7d74/tests/barcode

for f in cra_v2.bam cra_v2.bam.bai cra_v2_barcodes.txt; do
    curl -fsSL -o "$OUT/$f" "$BASE/$f"
done

printf 'sample,fastq_dir,bam,cell_barcodes\ncra_v2,,%s/cra_v2.bam,%s/cra_v2_barcodes.txt\n' "$OUT" "$OUT" > "$OUT/samplesheet.csv"
echo "Done: $OUT"; cat "$OUT/samplesheet.csv"
