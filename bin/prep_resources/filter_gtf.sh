#!/usr/bin/bash

# Usage: bash filter_gtf.sh input.gtf output.gtf

set -euo pipefail

if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <input.gtf> <output.gtf>"
  exit 1
fi

INPUT="$1"
OUTPUT="$2"

# Allowed biotypes (Cell Ranger–style whitelist)
ALLOWED_RE='^(protein_coding|lncRNA|antisense|IG_.*|TR_.*)$'

awk -v FS='\t' -v OFS='\t' -v allowed_re="$ALLOWED_RE" '
function get_attr(key, s,  pat, m) {
  pat = key " \"([^\"\n]+)\""
  if (match(s, pat, m)) return m[1]
  return ""
}
function get_biotype(attr,  bt) {
  bt = get_attr("gene_type", attr)
  if (bt == "") bt = get_attr("gene_biotype", attr)
  return bt
}

# ---------- Pass 1: collect allowed gene_ids ----------
FNR==NR {
  if ($0 ~ /^#/) next
  if ($3 == "gene") {
    bt  = get_biotype($9)
    gid = get_attr("gene_id", $9)
    if (gid != "" && bt ~ allowed_re) kept[gid]=1
  }
  next
}

# ---------- Pass 2: print header, allowed genes, and any feature of allowed genes ----------
{
  if ($0 ~ /^#/) { print; next }

  if ($3 == "gene") {
    bt  = get_biotype($9)
    if (bt ~ allowed_re) print
  } else {
    gid = get_attr("gene_id", $9)
    if (gid in kept) print
  }
}
' "$INPUT" "$INPUT" > "$OUTPUT"