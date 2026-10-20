#!/usr/bin/env python3

"""
Split a library BAM into batches of cells, in one pass.
Usage: split_bam_by_cb.py <in.bam> <barcodes.tsv> <batch_size|auto> [--barcode_tag CB]

Barcodes (first column of the file) are assigned to consecutive batches. Every read whose barcode tag is in the list
is written, in its original order, to the BAM of its batch; reads of other barcodes are dropped.
batch_size 'auto' targets about 75 batches per library, with 100 to 1000 cells each.
Writes batch_<NNNN>.bam and batch_<NNNN>.csv (the barcodes of the batch, one per line).
"""

import sys
import gzip
import math
import argparse
import pysam


##


parser = argparse.ArgumentParser(description='Split a BAM in batches of cells')
parser.add_argument('bam')
parser.add_argument('barcodes')
parser.add_argument('batch_size', help='cells per batch, or auto')
parser.add_argument('--barcode_tag', default='CB')
args = parser.parse_args()

MAX_BATCHES = 500       # open files at once
AUTO_TARGET, AUTO_MIN, AUTO_MAX = 75, 100, 1000


##


def batch_size_for(n_cells):
    if args.batch_size == 'auto':
        return min(max(math.ceil(n_cells / AUTO_TARGET), AUTO_MIN), AUTO_MAX)
    size = int(args.batch_size)
    if size < 1:
        sys.exit('split_bam_by_cb: batch_size must be >= 1')
    return size


def main():

    opener = gzip.open if args.barcodes.endswith('.gz') else open
    with opener(args.barcodes, 'rt') as fh:
        cells = list(dict.fromkeys(x.split('\t')[0].strip() for x in fh if x.strip()))
    if not cells:
        sys.exit(f'split_bam_by_cb: no barcodes in {args.barcodes}')

    size = batch_size_for(len(cells))
    n_batches = math.ceil(len(cells) / size)
    if n_batches > MAX_BATCHES:
        sys.exit(f'split_bam_by_cb: {n_batches} batches of {size} cells exceed the {MAX_BATCHES} allowed; raise --batch_size')
    batch_of = {c: i // size for i, c in enumerate(cells)}

    for b in range(n_batches):
        with open(f'batch_{b:04d}.csv', 'w') as fh:
            fh.writelines(f'{c}\n' for c in cells[b * size:(b + 1) * size])

    n_reads = [0] * n_batches
    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        outs = [pysam.AlignmentFile(f'batch_{b:04d}.bam', 'wb', template=bam) for b in range(n_batches)]
        try:
            for read in bam.fetch(until_eof=True):
                if not read.has_tag(args.barcode_tag):
                    continue
                b = batch_of.get(read.get_tag(args.barcode_tag))
                if b is not None:
                    outs[b].write(read)
                    n_reads[b] += 1
        finally:
            for o in outs:
                o.close()

    print(f'{len(cells)} barcodes -> {n_batches} batches of up to {size}; {sum(n_reads)} reads kept '
          f'(min {min(n_reads)}, max {max(n_reads)} per batch)')


##


if __name__ == '__main__':
    main()
