#!/usr/bin/python

"""
mgatk tenx counting for one chunk of cell barcodes of a library (mgatk 0.9.x).
Runs mgatk's own per-slice function (filter -> sort -> dedup -> per-cell, per-strand
allele counts), without its Snakemake wrapper.

Usage: mgatk_chunk.py <chunk.bam> <barcodes_<idx>.csv> <mt.fa> <out_prefix> [options]
Writes <out_prefix>.{A,C,G,T,coverage}.txt (+ <out_prefix>.depth.txt).
"""

import os
import argparse
from pathlib import Path
import pysam


##


parser = argparse.ArgumentParser(description='mgatk counting for one barcode chunk')
parser.add_argument('bam')
parser.add_argument('barcodes')
parser.add_argument('fasta', help='MT FASTA (a single contig)')
parser.add_argument('out_prefix')
parser.add_argument('--barcode_tag', default='CB')
parser.add_argument('--umi_tag', default='', help='UMI tag for duplicate removal (empty: none)')
parser.add_argument('--keep_duplicates', action='store_true')
parser.add_argument('--handle_overlap', action='store_true')
parser.add_argument('--proper_paired', action='store_true')
parser.add_argument('--base_qual', type=int, default=0)
parser.add_argument('--alignment_quality', type=int, default=0)
parser.add_argument('--NHmax', type=int, default=1)
parser.add_argument('--NMmax', type=int, default=4)
args = parser.parse_args()


##


def main():

    from mgatk.processing.one_slice_tenx import process_one_slice_tenx
    from mgatk.processing.barcodes import chunk_barcoded_bam

    # Reference
    with open(args.fasta) as fh:
        lines = [x.strip() for x in fh if x.strip()]
    names = [x[1:].split()[0] for x in lines if x.startswith('>')]
    if len(names) != 1:
        raise ValueError(f'{args.fasta}: expected a single contig, found {len(names)}')
    mito_chr = names[0]
    mito_length = sum(len(x) for x in lines if not x.startswith('>'))

    # mgatk expects its output tree
    out = Path('mgatk_out').resolve()
    for d in ['temp/barcode_files', 'temp/sparse_matrices', 'temp/ready_bam', 'temp/temp_bam',
              'logs/rmdupslogs', 'logs/filterlogs', 'qc/depth']:
        (out / d).mkdir(parents=True, exist_ok=True)
    sample = 'chunk'
    (out / 'temp/barcode_files' / f'{sample}.txt').write_text(Path(args.barcodes).read_text())

    config = {
        'output_directory': str(out),
        'mito_chr': mito_chr,
        'mito_length': mito_length,
        'fasta_file': args.fasta,
        'remove_duplicates': 'False' if args.keep_duplicates else 'True',
        'barcode_tag': args.barcode_tag,
        'umi_barcode': 'MU',
        'handle_overlap': 'True' if args.handle_overlap else 'False',
        'proper_paired': 'True' if args.proper_paired else 'False',
        'base_qual': args.base_qual,
        'alignment_quality': args.alignment_quality,
        'NHmax': args.NHmax,
        'NMmax': args.NMmax,
    }
    # As `mgatk tenx` does: write the chunk's reads with a synthetic "MU" tag (cell barcode + UMI), which is
    # what duplicate removal groups by. Without it, identical fragments of different cells would be merged.
    pysam.index(args.bam)
    barcoded = out / 'temp/barcoded_bams'
    barcoded.mkdir(parents=True, exist_ok=True)
    chunk_barcoded_bam(args.bam, str(barcoded), args.barcode_tag, args.barcodes, mito_chr, args.umi_tag or 'XX')
    chunk_bam = barcoded / (Path(args.barcodes).stem + '.bam')

    process_one_slice_tenx(config, str(chunk_bam), sample)

    for letter in ['A', 'C', 'G', 'T', 'coverage']:
        os.replace(out / f'temp/sparse_matrices/{sample}.{letter}.txt', f'{args.out_prefix}.{letter}.txt')
    os.replace(out / f'qc/depth/{sample}.depth.txt', f'{args.out_prefix}.depth.txt')


##


if __name__ == '__main__':
    main()
