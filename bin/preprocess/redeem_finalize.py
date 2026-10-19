#!/usr/bin/python

"""
Collapse the per-chunk RedeemV tables of a library into the tables MiTo reads, removing strand-biased variants.
Usage: redeem_finalize.py <out_dir> [--p_threshold 0.01] [--min_log2fc 1]
(chunk tables in cwd: *.RawGenotypes.<level>.txt and *.QualifiedTotalCts)

Strand bias (RedeemV's StrandBiases.R + RemoveStrandBias.py, with a plain binomial test):
 - per variant (position_ref_call), pooled over all the cells of the library, count the molecules
   seen on the plus strand and on the minus strand (columns 12 and 13 of the genotype rows);
 - two-sided binomial test of `plus` out of `plus + minus` against 0.5. The p-value is NOT corrected for the
   number of variants tested;
 - a variant is blacklisted if p < p_threshold and the strands differ by more than min_log2fc (log2 of the
   plus/minus ratio, 1 = 2-fold);
 - rows of blacklisted variants are removed, unless they have double-covered support (column 10, DB_Cts > 0).
Writes, in <out_dir>: RawGenotypes.<level>.StrandBalance.gz for the four levels, QualifiedTotalCts.gz,
StrandBiasTable.tsv.gz (every tested variant) and StrandBiasBlackList.
"""

import sys
import gzip
import glob
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest


##


parser = argparse.ArgumentParser()
parser.add_argument('out_dir')
parser.add_argument('--p_threshold', type=float, default=0.01)
parser.add_argument('--min_log2fc', type=float, default=1.0)
args = parser.parse_args()

LEVELS = ['Total', 'VerySensitive', 'Sensitive', 'Specific']
VARIANT, DB, PLUS, MINUS = 3, 9, 11, 12


##


def chunk_files(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f'No files match {pattern} in {Path.cwd()}')
    return files


def strand_table():
    """plus / minus molecule counts per variant over all cells (RawGenotypes.Total)."""
    plus, minus = {}, {}
    for path in chunk_files('*.RawGenotypes.Total.txt'):
        with open(path) as fh:
            for line in fh:
                f = line.split('\t')
                v = f[VARIANT]
                plus[v] = plus.get(v, 0) + int(float(f[PLUS]))
                minus[v] = minus.get(v, 0) + int(float(f[MINUS]))

    df = pd.DataFrame({'plus': pd.Series(plus, dtype='int64'), 'minus': pd.Series(minus, dtype='int64')})
    df.index.name = 'variant'
    df['N'] = df['plus'] + df['minus']
    df['p'] = [binomtest(int(k), int(n), 0.5).pvalue if n > 0 else 1.0 for k, n in zip(df['plus'], df['N'])]
    with np.errstate(divide='ignore', invalid='ignore'):
        df['abs_log2fc'] = np.abs(np.log2(df['plus'].to_numpy() / df['minus'].to_numpy()))   # inf if one strand is empty
    df['blacklisted'] = (df['p'] < args.p_threshold) & (df['abs_log2fc'] > args.min_log2fc)
    return df


##


def main():

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    table = strand_table()
    table.to_csv(out / 'StrandBiasTable.tsv.gz', sep='\t')
    blacklist = set(table.index[table['blacklisted']])
    with open(out / 'StrandBiasBlackList', 'w') as fh:
        fh.writelines(f'{v}\n' for v in sorted(blacklist))
    print(f'{len(table)} variants tested, {len(blacklist)} blacklisted (p<{args.p_threshold}, |log2fc|>{args.min_log2fc})')

    for level in LEVELS:
        n = removed = 0
        with gzip.open(out / f'RawGenotypes.{level}.StrandBalance.gz', 'wt', compresslevel=4) as o:
            for path in chunk_files(f'*.RawGenotypes.{level}.txt'):
                with open(path) as fh:
                    for line in fh:
                        n += 1
                        f = line.split('\t')
                        if f[VARIANT] in blacklist and float(f[DB]) == 0:
                            removed += 1
                            continue
                        o.write(line)
        print(f'{level}: {n} rows, {removed} removed')

    with gzip.open(out / 'QualifiedTotalCts.gz', 'wt', compresslevel=4) as o:
        for path in chunk_files('*.QualifiedTotalCts'):
            with open(path) as fh:
                for line in fh:
                    o.write(line)


##


if __name__ == '__main__':
    main()
