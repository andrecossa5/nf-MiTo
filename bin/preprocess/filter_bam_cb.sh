#!/usr/bin/env bash
# Keep the records of a BAM whose barcode tag is in a list of cell barcodes (one chunk of them).
# Usage: filter_bam_cb.sh <in.bam> <barcodes_<idx>.csv> [barcode_tag=CB]
# Writes filtered_<idx>.bam. Picard FilterSamReads; the barcodes go through an arguments file, so
# the number of barcodes is not limited by the length of the command line.
set -euo pipefail

bam=$1
cbs=$2
tag=${3:-CB}

[[ -s $bam ]] || { echo "filter_bam_cb: missing or empty BAM: $bam" >&2; exit 2; }
[[ -s $cbs ]] || { echo "filter_bam_cb: missing or empty barcode file: $cbs" >&2; exit 2; }

base=$(basename "$cbs")
idx=${base%%.*}
idx=${idx#*_}
out=filtered_${idx}.bam
args=filter_args_${idx}.txt

# One token per line: -TV <barcode>
awk 'NF { print "-TV"; print $1 }' "$cbs" > "$args"

picard FilterSamReads \
    -I "$bam" \
    -O "$out" \
    -TAG "$tag" \
    -FILTER includeTagValues \
    --arguments_file "$args" \
    --QUIET true \
    --COMPRESSION_LEVEL 1 \
    --MAX_RECORDS_IN_RAM 10000000

[[ -s $out ]] || { echo "filter_bam_cb: Picard produced no output ($out)" >&2; exit 3; }
rm -f "$args"
echo "$cbs: wrote $out ($(awk 'NF' "$cbs" | wc -l) barcodes, tag $tag)"
