#!/usr/bin/python

"""
RedeemV consensus calling for one chunk of cell barcodes (wraps RedeemV's mitoConsensus.py, unmodified).
Usage: redeem_consensus.py <chunk.bam> <barcodes_<idx>.csv> <mt.fa> <out_prefix> [options]

mitoConsensus.py works on a fixed folder layout, which is built here. Its output is completed as
RedeemV's Finalize does (AddQualifiedTotalCts): the depth at the cell and site, for the stringency level
of each file, is added as the 14th column of the genotype rows.

Writes <out_prefix>.RawGenotypes.<level>.txt (level: Total, VerySensitive, Sensitive, Specific)
and <out_prefix>.QualifiedTotalCts.
"""

import os
import sys
import shutil
import argparse
import subprocess
from pathlib import Path


##


parser = argparse.ArgumentParser(description='RedeemV mitoConsensus.py on a chunk')
parser.add_argument('bam')
parser.add_argument('barcodes')
parser.add_argument('fasta', help='MT FASTA (a single contig)')
parser.add_argument('out_prefix')
parser.add_argument('--redeemv_home', default=os.environ.get('REDEEMV_HOME', '/opt/redeemV'))
parser.add_argument('--barcode_tag', default='CB')
parser.add_argument('--base_quality', type=int, default=30)
args = parser.parse_args()

LEVELS = ['Total', 'VerySensitive', 'Sensitive', 'Specific']


##


def main():

    script = Path(args.redeemv_home) / 'mitoConsensus' / 'mitoConsensus.py'
    if not script.exists():
        raise FileNotFoundError(f'{script} not found (set --redeemv_home)')

    with open(args.fasta) as fh:
        seq = ''.join(x.strip() for x in fh if not x.startswith('>')).upper()

    sample = 'chunk'
    work = Path('redeem_work').resolve()
    for d in ['temp/barcode_files', 'temp/barcoded_bams', 'temp/sparse_matrices2.0', 'final']:
        (work / d).mkdir(parents=True, exist_ok=True)
    shutil.copy(args.barcodes, work / 'temp/barcode_files' / f'{sample}.txt')
    os.symlink(Path(args.bam).resolve(), work / 'temp/barcoded_bams' / f'{sample}.bam')
    with open(work / 'final/chrM_refAllele.txt', 'w') as fh:
        for pos, base in enumerate(seq, start=1):
            fh.write(f'{pos}\t{base}\n')

    subprocess.run([sys.executable, str(script), sample, str(work), args.barcode_tag, str(args.base_quality)], check=True)

    out = work / 'temp/sparse_matrices2.0'
    # Depth at (cell, site), per level: QualifiedTotalCts is  cell, pos, Total, VerySensitive, Sensitive, Specific
    depth = {}
    with open(out / f'{sample}.QualifiedTotalCts') as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            depth[(f[0], f[1])] = f[2:6]
    shutil.copy(out / f'{sample}.QualifiedTotalCts', f'{args.out_prefix}.QualifiedTotalCts')

    for k, level in enumerate(LEVELS):
        with open(out / f'{sample}.RawGenotypes.{level}') as fh, open(f'{args.out_prefix}.RawGenotypes.{level}.txt', 'w') as o:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                o.write('\t'.join(f + [depth[(f[1], f[2])][k]]) + '\n')

    shutil.rmtree(work)


##


if __name__ == '__main__':
    main()
