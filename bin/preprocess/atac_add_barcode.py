#!/usr/bin/python

"""
Attach the cell barcode of a 10x scATAC library to the genomic reads.
Usage: atac_add_barcode.py <R1.fastq.gz> <R2_barcode.fastq.gz> <R3.fastq.gz> <cell_barcodes.txt> <out_prefix>

The barcode read (R2, 16 nt) is matched to the list of valid cell barcodes (one per line, with or
without a "-1" suffix): exact match, or a single unambiguous mismatch. Read pairs whose barcode
does not match a listed cell are dropped (not needed downstream, and this makes mapping much
cheaper). Writes <out_prefix>_R1.fastq.gz and <out_prefix>_R3.fastq.gz with the corrected barcode
appended to the read name as "<name>_BC:<barcode>" (the barcode as listed, suffix included).
"""

import sys
import gzip
import itertools


##


path_r1, path_r2, path_r3, path_cbs, out_prefix = sys.argv[1:6]
BC_LEN = 16


##


def neighbours(bc):
    """All sequences at Hamming distance 1 (substitutions of A,C,G,T,N)."""
    for i in range(len(bc)):
        for b in 'ACGTN':
            if b != bc[i]:
                yield bc[:i] + b + bc[i+1:]


def read_barcodes(path):
    listed = {}
    with open(path) as fh:
        for line in fh:
            bc = line.strip()
            if bc:
                listed[bc[:BC_LEN]] = bc
    if not listed:
        raise ValueError(f'No barcodes in {path}')
    return listed


def build_lookup(listed):
    """Barcode read -> listed barcode: exact matches, plus unambiguous 1-mismatch ones."""
    lookup = dict(listed)
    seen = {}
    for core in listed:
        for nb in neighbours(core):
            if nb in listed:
                continue
            seen[nb] = None if nb in seen else core   # ambiguous neighbours are discarded
    lookup.update({nb: listed[core] for nb, core in seen.items() if core is not None})
    return lookup


def clean_name(head):
    """'@name/1 comment' -> '@name' (the mates of a pair must share the name)."""
    name = head.split()[0]
    return name[:-2] if len(name) > 3 and name[-2] == '/' and name[-1].isdigit() else name


def records(handle):
    while True:
        head = handle.readline()
        if not head:
            return
        yield head, handle.readline(), handle.readline(), handle.readline()


##


def main():

    lookup = build_lookup(read_barcodes(path_cbs))
    n = kept = exact = 0

    with gzip.open(path_r1, 'rt') as f1, gzip.open(path_r2, 'rt') as f2, gzip.open(path_r3, 'rt') as f3, \
         gzip.open(f'{out_prefix}_R1.fastq.gz', 'wt', compresslevel=3) as o1, \
         gzip.open(f'{out_prefix}_R3.fastq.gz', 'wt', compresslevel=3) as o3:

        for a, b, c in itertools.zip_longest(records(f1), records(f2), records(f3)):
            if a is None or b is None or c is None:
                raise ValueError('R1, R2 and R3 have different numbers of reads')
            n += 1
            core = b[1].strip()[:BC_LEN]
            bc = lookup.get(core)
            if bc is None:
                continue
            kept += 1
            exact += core == bc[:BC_LEN]
            name = clean_name(a[0])
            o1.write(f'{name}_BC:{bc}\n{a[1]}+\n{a[3]}')
            o3.write(f'{name}_BC:{bc}\n{c[1]}+\n{c[3]}')

    print(f'{n} read pairs, {kept} with a listed cell barcode ({exact} exact)')


##


if __name__ == '__main__':
    main()
