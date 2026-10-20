#!/usr/bin/env python3

"""
mgatk `bcall` counting for a batch of cells, one cell at a time.
Each cell gets its own BAM (all its reads, MT contig); mgatk's per-sample function (filter -> sort -> dedup ->
per-strand allele counts -> depth) runs on it, without the Snakemake wrapper. mgatk is imported once per batch.

Usage: mgatk_batch.py <batch.bam> <batch_id> <mt.fa> [options]
Writes <batch_id>.{A,C,G,T,coverage,depth}.txt: the tables of the cells of the batch, concatenated.
"""

import os
import shutil
import argparse
import subprocess
from pathlib import Path


##


parser = argparse.ArgumentParser(description='mgatk bcall counting for a batch of cells')
parser.add_argument('bam')
parser.add_argument('batch_id')
parser.add_argument('fasta', help='MT FASTA (a single contig)')
parser.add_argument('--barcode_tag', default='CB')
parser.add_argument('--umi_tag', default='', help='UMI tag for duplicate removal (empty: coordinates only)')
parser.add_argument('--keep_duplicates', action='store_true')
parser.add_argument('--handle_overlap', action='store_true')
parser.add_argument('--proper_paired', action='store_true')
parser.add_argument('--base_qual', type=int, default=0)
parser.add_argument('--alignment_quality', type=int, default=0)
parser.add_argument('--NHmax', type=int, default=1)
parser.add_argument('--NMmax', type=int, default=4)
args = parser.parse_args()

TABLES = ['A', 'C', 'G', 'T', 'coverage', 'depth']


##


def main():

    from mgatk.processing.one_sample import process_one_sample

    with open(args.fasta) as fh:
        lines = [x.strip() for x in fh if x.strip()]
    names = [x[1:].split()[0] for x in lines if x.startswith('>')]
    if len(names) != 1:
        raise ValueError(f'{args.fasta}: expected a single contig, found {len(names)}')
    mito_length = sum(len(x) for x in lines if not x.startswith('>'))

    cells_dir = Path('cells')
    cells_dir.mkdir()
    subprocess.run(['samtools', 'split', '-M', '-1', '-d', args.barcode_tag, '-f', 'cells/%!.bam', args.bam], check=True)
    for stray in cells_dir.glob('-.bam'):               # reads without the barcode tag
        stray.unlink()

    out = Path('mgatk_out').resolve()
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

    handles = {t: open(f'{args.batch_id}.{t}.txt', 'wb') for t in TABLES}
    cell_bams = sorted(cells_dir.glob('*.bam'))
    for k, cell_bam in enumerate(cell_bams, start=1):
        cell = cell_bam.name[:-len('.bam')]
        for d in ['temp/ready_bam', 'temp/temp_bam', 'temp/sparse_matrices', 'logs/rmdupslogs', 'logs/filterlogs', 'qc/depth']:
            (out / d).mkdir(parents=True, exist_ok=True)
        process_one_sample(config, str(cell_bam.resolve()), str(out / 'temp/ready_bam' / f'{cell}.qc.bam'), cell)
        for t in TABLES:
            src = out / (f'qc/depth/{cell}.depth.txt' if t == 'depth' else f'temp/sparse_matrices/{cell}.{t}.txt')
            with open(src, 'rb') as fh:
                shutil.copyfileobj(fh, handles[t])
        shutil.rmtree(out)                              # per-cell intermediates
        cell_bam.unlink()
        if k % 50 == 0:
            print(f'{k}/{len(cell_bams)} cells', flush=True)
    for h in handles.values():
        h.close()
    print(f'{args.batch_id}: {len(cell_bams)} cells')


##


if __name__ == '__main__':
    main()
