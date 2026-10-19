#!/usr/bin/python

"""
Turn the bwa-mem2 SAM stream of barcoded scATAC reads into the MT BAM MiTo's tools need.
Usage: bwa-mem2 mem ... | atac_tag_bam.py <mt_contig> <min_mapq> <barcode_tag> <out.bam>

Keeps properly paired, primary alignments on the MT contig with MAPQ >= min_mapq, moves the barcode
that atac_add_barcode.py appended to the read name ("_BC:<barcode>") into <barcode_tag>, and writes
an unsorted BAM (sort and index it afterwards).
"""

import sys
import pysam


##


mt_contig, min_mapq, barcode_tag, out_path = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]


##


def main():

    n = kept = 0
    with pysam.AlignmentFile('-', 'r') as sam:
        header = sam.header.to_dict()
        # keep the whole reference in the header (a later sort/index and MiTo need it), but only MT reads
        with pysam.AlignmentFile(out_path, 'wb', header=header) as out:
            for read in sam.fetch(until_eof=True):
                n += 1
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                if read.reference_name != mt_contig or read.mapping_quality < min_mapq or not read.is_proper_pair:
                    continue
                name, _, bc = read.query_name.rpartition('_BC:')
                if not name:
                    continue
                read.query_name = name
                read.set_tag(barcode_tag, bc, value_type='Z')
                out.write(read)
                kept += 1

    print(f'{n} alignments, {kept} kept (MT, proper pair, MAPQ>={min_mapq})', file=sys.stderr)


##


if __name__ == '__main__':
    main()
