#!/usr/bin/python

"""
Translate Cell Ranger ARC (RNA-space) cell barcodes to the ATAC-space barcodes found in the reads
of the ATAC library of a 10x Multiome experiment.
Usage: rna2atac_barcodes.py <rna_cells.tsv(.gz)> <WhiteList.RNA> <WhiteList.ATAC> <out.tsv>

In a Multiome bead the RNA and the ATAC barcode are linked: the two whitelists are line-aligned (line i of
one pairs with line i of the other). Cell Ranger ARC reports cells with the RNA barcode ("<16 nt>-1"), while
the ATAC reads carry the ATAC barcode. Writes "<RNA barcode as listed>\t<ATAC barcode>" per cell, so that
the cells keep the names of the RNA side (and of the rest of the Multiome object).
"""

import sys
import gzip


##


rna_cells, wl_rna, wl_atac, out = sys.argv[1:5]


##


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def main():

    with opener(wl_rna) as r, opener(wl_atac) as a:
        rna = [x.strip() for x in r if x.strip()]
        atac = [x.strip() for x in a if x.strip()]
    if len(rna) != len(atac):
        raise ValueError(f'The whitelists have different sizes: {len(rna)} vs {len(atac)}')
    translate = dict(zip(rna, atac))

    n = missing = 0
    with opener(rna_cells) as fh, open(out, 'w') as o:
        for line in fh:
            name = line.strip().split('\t')[0]
            if not name:
                continue
            n += 1
            bc = translate.get(name[:16])
            if bc is None:
                missing += 1
                continue
            o.write(f'{name}\t{bc}\n')

    print(f'{n} barcodes, {n - missing} translated, {missing} not in the whitelist', file=sys.stderr)
    if n == missing:
        raise ValueError('No barcode found in the RNA whitelist: are these RNA-space barcodes?')


##


if __name__ == '__main__':
    main()
