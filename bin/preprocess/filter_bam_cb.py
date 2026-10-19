#!/usr/bin/python

"""
Keep the records of a BAM whose barcode tag is in a list of cell barcodes (one chunk of them).
Usage: filter_bam_cb.py <in.bam> <barcodes_<idx>.csv> [barcode_tag=CB]
Writes filtered_<idx>.bam, preserving the record order (and so the sorting) of the input.
"""

import sys
import pysam


##


i = sys.argv[1]
path_cbs = sys.argv[2]
tag = sys.argv[3] if len(sys.argv) > 3 else 'CB'


##


def main():

    idx = path_cbs.split('.')[0].split('_')[1]
    with open(path_cbs) as fh:
        cbs = {x.strip() for x in fh if x.strip()}

    n_in = n_out = 0
    with pysam.AlignmentFile(i, 'rb') as bam_in, \
         pysam.AlignmentFile(f'filtered_{idx}.bam', 'wb', template=bam_in) as bam_out:
        for read in bam_in.fetch(until_eof=True):
            n_in += 1
            if read.has_tag(tag) and read.get_tag(tag) in cbs:
                bam_out.write(read)
                n_out += 1

    print(f'{path_cbs}: kept {n_out}/{n_in} records ({len(cbs)} barcodes, tag {tag})')


##


if __name__ == '__main__':
    main()
