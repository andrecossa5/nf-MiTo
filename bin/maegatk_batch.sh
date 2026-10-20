#!/usr/bin/env bash
# maegatk for a batch of cells: split the batch BAM by cell, run oneSample_maegatk.py on each cell, concatenate the tables.
# Usage: maegatk_batch.sh <batch.bam> <batch_id> <mt.fa> <barcode_tag> [oneSample_maegatk.py options]
# Writes <batch_id>.{A,C,G,T,coverage,depth}.txt
set -euo pipefail

bam=$1
batch=$2
fasta=$(realpath "$3")
tag=$4
shift 4

[[ -s $bam ]] || { echo "maegatk_batch: missing or empty BAM: $bam" >&2; exit 2; }

mkdir -p cells work
samtools split -M -1 -d "$tag" -f 'cells/%!.bam' "$bam"
rm -f -- cells/-.bam                                   # reads without the barcode tag

tables="A C G T coverage depth"
for t in $tables; do : > "${batch}.${t}.txt"; done

n=0
total=$(ls cells | wc -l)
for cell_bam in cells/*.bam; do
    [[ -e $cell_bam ]] || break
    cell=$(basename "$cell_bam" .bam)
    mkdir -p "work/$cell"
    ( cd "work/$cell" && oneSample_maegatk.py "$(realpath "../../$cell_bam")" "$cell" "$fasta" "$@" )
    for t in $tables; do cat "work/$cell/$cell.$t.txt" >> "${batch}.${t}.txt"; done
    rm -rf "work/$cell" "$cell_bam"
    n=$((n + 1))
    (( n % 50 == 0 )) && echo "$n/$total cells" || true
done
echo "$batch: $n cells"
