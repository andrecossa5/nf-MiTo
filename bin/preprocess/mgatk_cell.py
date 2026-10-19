#!/usr/bin/python

"""
mgatk `bcall` counting for one cell: its own BAM (all reads of the cell, MT contig).
Runs mgatk's per-sample function (filter -> sort -> dedup -> per-strand allele counts -> depth),
without its Snakemake wrapper.

Usage: mgatk_cell.py <cell.bam> <cell> <mt.fa> [options]
Writes <cell>.{A,C,G,T,coverage,depth}.txt
"""

import os
import argparse
from pathlib import Path


##


parser = argparse.ArgumentParser(description='mgatk bcall counting for one cell')
parser.add_argument('bam')
parser.add_argument('cell')
parser.add_argument('fasta', help='MT FASTA (a single contig)')
parser.add_argument('--umi_tag', default='', help='UMI tag for duplicate removal (empty: coordinates only)')
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

    from mgatk.processing.one_sample import process_one_sample

    with open(args.fasta) as fh:
        lines = [x.strip() for x in fh if x.strip()]
    names = [x[1:].split()[0] for x in lines if x.startswith('>')]
    if len(names) != 1:
        raise ValueError(f'{args.fasta}: expected a single contig, found {len(names)}')
    mito_length = sum(len(x) for x in lines if not x.startswith('>'))

    out = Path('mgatk_out').resolve()
    for d in ['temp/ready_bam', 'temp/temp_bam', 'temp/sparse_matrices', 'logs/rmdupslogs', 'logs/filterlogs', 'qc/depth']:
        (out / d).mkdir(parents=True, exist_ok=True)

    config = {
        'output_directory': str(out),
        'mito_chr': names[0],
        'mito_length': mito_length,
        'fasta_file': args.fasta,
        'remove_duplicates': 'False' if args.keep_duplicates else 'True',
        'umi_barcode': args.umi_tag,
        'emit_base_qualities': False,
        'handle_overlap': 'True' if args.handle_overlap else 'False',
        'proper_paired': 'True' if args.proper_paired else 'False',
        'base_qual': args.base_qual,
        'alignment_quality': args.alignment_quality,
        'NHmax': args.NHmax,
        'NMmax': args.NMmax,
    }
    output_bam = str(out / 'temp/ready_bam' / f'{args.cell}.qc.bam')
    process_one_sample(config, str(Path(args.bam).resolve()), output_bam, args.cell)

    for letter in ['A', 'C', 'G', 'T', 'coverage']:
        os.replace(out / f'temp/sparse_matrices/{args.cell}.{letter}.txt', f'{args.cell}.{letter}.txt')
    os.replace(out / f'qc/depth/{args.cell}.depth.txt', f'{args.cell}.depth.txt')


##


if __name__ == '__main__':
    main()
