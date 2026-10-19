#!/usr/bin/python

"""
Attach the cell barcode of a 10x ATAC / Multiome-ATAC library to the genomic reads.
Usage: atac_add_barcode.py <R1.fastq.gz> <R2_barcode.fastq.gz> <R3.fastq.gz> <cell_barcodes> <out_prefix> [--rule auto]

<cell_barcodes>: one barcode per line (with or without a "-1" suffix) or two tab-separated columns,
"<name>\t<barcode as in the reads>". With two columns, the second is matched to the reads and the first is
what is written (e.g. RNA-space names for an ATAC library, see rna2atac_barcodes.py).

Where the 16-nt barcode sits in the barcode read (offset, orientation) depends on the kit and the
sequencer (16 nt: 10x scATAC; 24 nt: Multiome i5, usually reverse complement of bases 8-24). With
--rule auto it is chosen by sampling reads against the cell list; --rule fwd:<offset> / rc:<offset> forces it.

The barcode is matched to the list exactly or with one unambiguous mismatch. Read pairs whose barcode is not
a listed cell are dropped (not needed downstream, and mapping gets much cheaper). Writes
<out_prefix>_R1.fastq.gz and <out_prefix>_R3.fastq.gz with "_BC:<name>" appended to the read name.
"""

import sys
import gzip
import argparse
import itertools


##


parser = argparse.ArgumentParser()
parser.add_argument('r1')
parser.add_argument('r2')
parser.add_argument('r3')
parser.add_argument('cell_barcodes')
parser.add_argument('out_prefix')
parser.add_argument('--rule', default='auto', help='auto | fwd:<offset> | rc:<offset>')
parser.add_argument('--detect_reads', type=int, default=50000)
parser.add_argument('--min_match_rate', type=float, default=0.02)
args = parser.parse_args()

BC_LEN = 16
COMP = str.maketrans('ACGTN', 'TGCAN')


##


def revcomp(s):
    return s[::-1].translate(COMP)


def neighbours(bc):
    """All sequences at Hamming distance 1 (substitutions of A,C,G,T,N)."""
    for i in range(len(bc)):
        for b in 'ACGTN':
            if b != bc[i]:
                yield bc[:i] + b + bc[i+1:]


def read_listed(path):
    """{barcode as in the reads (16 nt): name written in the read names}"""
    listed = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not parts[0].strip():
                continue
            name = parts[0].strip()
            key = (parts[1].strip() if len(parts) > 1 else name)[:BC_LEN]
            listed[key] = name
    if not listed:
        raise ValueError(f'No barcodes in {path}')
    return listed


def build_lookup(listed):
    """Barcode read -> name: exact matches, plus unambiguous 1-mismatch ones."""
    lookup = dict(listed)
    seen = {}
    for key in listed:
        for nb in neighbours(key):
            if nb in listed:
                continue
            seen[nb] = None if nb in seen else key   # ambiguous neighbours are discarded
    lookup.update({nb: listed[key] for nb, key in seen.items() if key is not None})
    return lookup


def window(seq, rule):
    mode, offset = rule
    w = seq[offset:offset + BC_LEN]
    return revcomp(w) if mode == 'rc' else w


def records(handle):
    while True:
        head = handle.readline()
        if not head:
            return
        yield head, handle.readline(), handle.readline(), handle.readline()


def clean_name(head):
    """'@name/1 comment' -> '@name' (the mates of a pair must share the name)."""
    name = head.split()[0]
    return name[:-2] if len(name) > 3 and name[-2] == '/' and name[-1].isdigit() else name


def choose_rule(listed):
    """Barcode window (orientation, offset) that best matches the cell list in the first reads."""
    if args.rule != 'auto':
        mode, offset = args.rule.split(':')
        return mode, int(offset)

    seqs = []
    with gzip.open(args.r2, 'rt') as fh:
        for rec in itertools.islice(records(fh), args.detect_reads):
            seqs.append(rec[1].strip())
    if not seqs:
        raise ValueError(f'{args.r2} is empty')
    length = len(seqs[0])
    if length < BC_LEN:
        raise ValueError(f'Barcode read of {length} nt, shorter than {BC_LEN}')

    rates = {}
    for mode in ['fwd', 'rc']:
        for offset in range(length - BC_LEN + 1):
            hits = sum(window(s, (mode, offset)) in listed for s in seqs)
            rates[(mode, offset)] = hits / len(seqs)
    best = max(rates, key=rates.get)
    if rates[best] < args.min_match_rate:
        top = sorted(rates.items(), key=lambda x: -x[1])[:3]
        raise ValueError(f'No barcode window matches the cell list (best: {top}). Check the list and its barcode space')
    print(f'Barcode window: {best[0]}:{best[1]} ({rates[best]:.1%} of {len(seqs)} reads match a listed cell)', file=sys.stderr)
    return best


##


def main():

    listed = read_listed(args.cell_barcodes)
    rule = choose_rule(listed)
    lookup = build_lookup(listed)
    n = kept = exact = 0

    with gzip.open(args.r1, 'rt') as f1, gzip.open(args.r2, 'rt') as f2, gzip.open(args.r3, 'rt') as f3, \
         gzip.open(f'{args.out_prefix}_R1.fastq.gz', 'wt', compresslevel=3) as o1, \
         gzip.open(f'{args.out_prefix}_R3.fastq.gz', 'wt', compresslevel=3) as o3:

        for a, b, c in itertools.zip_longest(records(f1), records(f2), records(f3)):
            if a is None or b is None or c is None:
                raise ValueError('R1, R2 and R3 have different numbers of reads')
            n += 1
            key = window(b[1].strip(), rule)
            name_bc = lookup.get(key)
            if name_bc is None:
                continue
            kept += 1
            exact += key in listed
            name = clean_name(a[0])
            o1.write(f'{name}_BC:{name_bc}\n{a[1]}+\n{a[3]}')
            o3.write(f'{name}_BC:{name_bc}\n{c[1]}+\n{c[3]}')

    print(f'{n} read pairs, {kept} with a listed cell barcode ({exact} exact)')


##


if __name__ == '__main__':
    main()
