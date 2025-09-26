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

# Simple approach: keep protein_coding genes and basic gene types
# First pass: copy header and identify allowed genes
grep '^#' "$INPUT" > "$OUTPUT" || true

# Second pass: filter for allowed biotypes
awk -F'\t' '
BEGIN { OFS="\t" }
!/^#/ {
    # Extract gene_type or gene_biotype from attributes column
    if (match($9, /gene_type "([^"]+)"/, arr)) {
        biotype = arr[1]
    } else if (match($9, /gene_biotype "([^"]+)"/, arr)) {
        biotype = arr[1]
    } else {
        biotype = ""
    }
    
    # Keep protein coding, lncRNA, antisense, and immunoglobulin genes
    if (biotype ~ /^(protein_coding|lncRNA|antisense|IG_.*|TR_.*)$/) {
        print $0
    }
}' "$INPUT" >> "$OUTPUT"

echo "GTF filtering completed"