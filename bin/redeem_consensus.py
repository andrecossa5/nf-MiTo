#!/usr/bin/env python3

"""
RedeemV molecule consensus for a batch of cells, streamed one cell at a time.
Usage: redeem_consensus.py <batch.bam> <barcodes.csv> <mt.fa> <out_prefix> [options]

Same rules, and same output, as RedeemV's mitoConsensus.py (chenweng1991/redeemV @ 5c34940) followed by its
AddQualifiedTotalCts step, with a memory use set by the largest cell instead of by the whole BAM:
 - the BAM is sorted by the barcode tag (samtools), so the reads of a cell are contiguous;
 - a molecule is (cell, start, end) of the forward read of a read pair; pairs are formed by read name inside the cell
   (names with a number of records other than 2, and pairs with both reads on the same strand, are discarded);
 - per molecule, per position: single-strand (SG) and double-strand (DB, mates overlapping and agreeing) base counts,
   bases kept if quality > base_quality, strand flags; the base with most counts is the call;
 - the four stringency levels (Total, VerySensitive, Sensitive, Specific) from family size and consensus fraction.
Genotypes are held as sparse dicts (covered positions only), not as 16569 x 5 arrays.

Writes <out_prefix>.RawGenotypes.<level>.txt (13 columns of mitoConsensus.py + the depth at the cell and site for
that level as 14th) and <out_prefix>.QualifiedTotalCts (cell, pos, Total, VerySensitive, Sensitive, Specific; all
positions of all cells of the barcode file).
"""

import os
import argparse
import subprocess
import tempfile
from itertools import groupby
import numpy as np
import pysam


##


parser = argparse.ArgumentParser(description='RedeemV molecule consensus, streamed by cell')
parser.add_argument('bam')
parser.add_argument('barcodes')
parser.add_argument('fasta', help='MT FASTA (a single contig)')
parser.add_argument('out_prefix')
parser.add_argument('--barcode_tag', default='CB')
parser.add_argument('--base_quality', type=int, default=30)
parser.add_argument('--threads', type=int, default=1)
args = parser.parse_args()

LEVELS = ['Total', 'VerySensitive', 'Sensitive', 'Specific']
LETTERS = 'ACGTN'
IDX = {b: i for i, b in enumerate(LETTERS)}


##


def read_reference(path):
    with open(path) as fh:
        return ''.join(x.strip() for x in fh if not x.startswith('>')).upper()


def sort_by_cell(bam, tag, threads, tmpdir):
    out = os.path.join(tmpdir, 'by_cell.bam')
    subprocess.run(['samtools', 'sort', '-t', tag, '-@', str(threads), '-m', '1G',
                    '-T', os.path.join(tmpdir, 'sort'), '-o', out, bam], check=True)
    return out


def molecules(reads, tag):
    """Molecule -> read-pair names, and name -> records, for the reads of one cell (as mitoConsensus.py)."""
    pairs = {}
    for read in reads:
        pairs.setdefault(read.query_name, []).append(read)
    mols = {}
    for name, rr in pairs.items():
        if len(rr) != 2:
            continue
        r0, r1 = rr
        if r0.is_reverse and not r1.is_reverse:
            fwd = r1
        elif not r0.is_reverse and r1.is_reverse:
            fwd = r0
        else:
            continue
        if fwd.has_tag(tag):
            cell = fwd.get_tag(tag)
            if fwd.template_length > 0:
                key = (cell, fwd.reference_start, fwd.reference_start + fwd.template_length)
            else:
                key = (cell, fwd.reference_start - fwd.template_length, fwd.reference_start)
            mols.setdefault(key, []).append(rr)
    return mols


def genotype_molecule(read_pairs, min_q):
    """Single-strand / double-strand base counts and strand flags per position of one molecule (sparse)."""
    sg, db, strand = {}, {}, {}

    def add(table, pos, k):
        v = table.get(pos)
        if v is None:
            v = table[pos] = [0, 0, 0, 0, 0]
        v[k] += 1

    def add_strand(pos, k):
        v = strand.get(pos)
        if v is None:
            v = strand[pos] = [0, 0]
        v[k] += 1

    for r0, r1 in read_pairs:
        p0 = r0.get_aligned_pairs(matches_only=True)
        p1 = r1.get_aligned_pairs(matches_only=True)
        s0, s1 = r0.query_sequence, r1.query_sequence
        q0, q1 = r0.query_qualities, r1.query_qualities
        overlap = {r for _, r in p0} & {r for _, r in p1}
        for q, r in p0:
            if r in overlap:
                continue
            if q0[q] > min_q:
                add(sg, r, IDX[s0[q]])
                add_strand(r, int(r0.is_reverse))
            else:
                add(sg, r, 4)
        for (qa, ra), (qb, rb) in zip([x for x in p0 if x[1] in overlap], [x for x in p1 if x[1] in overlap]):
            if s0[qa] == s1[qb]:
                if q0[qa] > min_q or q1[qb] > min_q:
                    add(db, ra, IDX[s0[qa]])
                    add_strand(ra, 0)
                    add_strand(ra, 1)
                else:
                    add(db, ra, 4)
            else:
                add(db, ra, 4)
        for q, r in p1:
            if r in overlap:
                continue
            if q1[q] > min_q:
                add(sg, r, IDX[s1[q]])
                add_strand(r, int(r1.is_reverse))
            else:
                add(sg, r, 4)
    return sg, db, strand


def f(x):
    return str(float(x))


def call_cell(cell, reads, seq, min_q, tag, rows):
    """Genotype rows (per level) and depth counts (n_levels x positions) of one cell."""
    counts = np.zeros((4, len(seq)), dtype=np.int64)       # molecules per level and position
    mols = molecules(reads, tag)
    zero = [0, 0, 0, 0, 0]
    for key, read_pairs in mols.items():
        m = f'{key[0]}_{key[1]}_{key[2]}'
        sg, db, strand = genotype_molecule(read_pairs, min_q)
        for i in sorted(set(sg) | set(db)):
            a, b = sg.get(i, zero), db.get(i, zero)
            cur = [a[k] + b[k] for k in range(4)]
            fam = sum(cur)
            if fam == 0:
                continue
            ci = cur.index(max(cur))
            call = LETTERS[ci]
            ref = seq[i]
            gt, sgc, dbc = cur[ci], a[ci], b[ci]
            css = gt / fam
            st = strand.get(i, (0, 0))
            counts[0, i] += 1
            out = '\t'.join([m, cell, str(i + 1), f'{i + 1}_{ref}_{call}', call, ref, f(fam), f(gt), str(css),
                             f(dbc), f(sgc), str(int(st[0] > 0)), str(int(st[1] > 0))])
            if call != ref:
                rows[0].append(out)
            if dbc == 0:
                thr = [(css > 0.75 and fam >= 2), (css > 0.75 and fam >= 3), (css > 0.9 and fam >= 4)]
            else:
                thr = [(css > 0.75 and fam >= 1), (css > 0.75 and fam >= 2), (css > 0.9 and fam >= 3)]
            for k, ok in enumerate(thr, start=1):
                if ok:
                    counts[k, i] += 1
                    if call != ref:
                        rows[k].append(out)
    return counts


def write_depth(fh, cell, counts):
    nz = {int(i) for i in np.flatnonzero(counts.any(axis=0))}
    zero = '0.0\t0.0\t0.0\t0.0'
    fh.writelines(
        f'{cell}\t{i + 1}\t' + (f'{counts[0, i]}.0\t{counts[1, i]}.0\t{counts[2, i]}.0\t{counts[3, i]}.0' if i in nz else zero) + '\n'
        for i in range(counts.shape[1])
    )


##


def main():

    seq = read_reference(args.fasta)
    with open(args.barcodes) as fh:
        barcodes = [x.strip() for x in fh if x.strip()]
    seen = set()

    geno = [open(f'{args.out_prefix}.RawGenotypes.{l}.txt', 'w') for l in LEVELS]
    depth = open(f'{args.out_prefix}.QualifiedTotalCts', 'w')

    with tempfile.TemporaryDirectory(dir='.') as tmp:
        bam = sort_by_cell(args.bam, args.barcode_tag, args.threads, tmp)
        with pysam.AlignmentFile(bam, 'rb') as af:
            tagged = (r for r in af if r.has_tag(args.barcode_tag))
            for cell, reads in groupby(tagged, key=lambda r: r.get_tag(args.barcode_tag)):
                rows = [[], [], [], []]
                counts = call_cell(cell, reads, seq, args.base_quality, args.barcode_tag, rows)
                seen.add(cell)
                for k in range(4):                      # depth at the cell and site for the level, 14th column
                    for line in rows[k]:
                        pos = int(line.split('\t', 3)[2]) - 1
                        geno[k].write(f'{line}\t{counts[k, pos]}.0\n')
                write_depth(depth, cell, counts)
        for cell in barcodes:                           # cells without molecules: zero depth
            if cell not in seen:
                write_depth(depth, cell, np.zeros((4, len(seq)), dtype=np.int64))

    for fh in geno + [depth]:
        fh.close()


##


if __name__ == '__main__':
    main()
