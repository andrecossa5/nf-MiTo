#!/usr/bin/env python3

"""
Attach the cell barcode of a 10x ATAC / Multiome-ATAC library to the genomic reads.
Usage: atac_add_barcode.py <R1.fastq.gz> <R2_barcode.fastq.gz> <R3.fastq.gz> <cell_barcodes> <out_prefix> [--rule auto] [--compress_threads 4]

<cell_barcodes>: one barcode per line (with or without a "-1" suffix) or two tab-separated columns,
"<name>\t<barcode as in the reads>". With two columns, the second is matched to the reads and the first is
what is written (e.g. RNA-space names for an ATAC library, see rna2atac_barcodes.py).

Where the 16-nt barcode sits in the barcode read (offset, orientation) depends on the kit and the
sequencer (16 nt: 10x scATAC; 24 nt: Multiome i5, usually reverse complement of bases 8-24). With
--rule auto it is chosen by sampling reads against the cell list; --rule fwd:<offset> / rc:<offset> forces it.

The barcode is matched to the list exactly or with one unambiguous mismatch. Read pairs whose barcode is not
a listed cell are dropped (not needed downstream, and mapping gets much cheaper). Writes
<out_prefix>_R1.fastq.gz and <out_prefix>_R3.fastq.gz with "_BC:<name>" appended to the read name.

Speed: reads are handled as bytes, and gzip is done by external pigz processes (gzip if pigz is not installed): the
three inputs are decompressed and the two outputs compressed (--compress_threads each) in parallel with this process.
"""

import sys
import shutil
import argparse
import itertools
import subprocess
import gzip


##


parser = argparse.ArgumentParser()
parser.add_argument('r1')
parser.add_argument('r2')
parser.add_argument('r3')
parser.add_argument('cell_barcodes')
parser.add_argument('out_prefix')
parser.add_argument('--rule', default='auto', help='auto | fwd:<offset> | rc:<offset>')
parser.add_argument('--detect_reads', type=int, default=50000)
parser.add_argument('--min_match_rate', type=float, default=0.02)
parser.add_argument('--compress_threads', type=int, default=4, help='pigz threads for each of the two outputs')
args = parser.parse_args()

BC_LEN = 16
COMP = bytes.maketrans(b'ACGTN', b'TGCAN')
BUFSIZE = 1 << 20
PIGZ = shutil.which('pigz')


##


def revcomp(s):
    return s[::-1].translate(COMP)


def neighbours(bc):
    """All sequences at Hamming distance 1 (substitutions of A,C,G,T,N)."""
    for i in range(len(bc)):
        for b in b'ACGTN':
            if b != bc[i]:
                yield bc[:i] + bytes((b,)) + bc[i+1:]


def read_listed(path):
    """{barcode as in the reads (16 nt): name written in the read names}"""
    listed = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not parts[0].strip():
                continue
            name = parts[0].strip()
            key = (parts[1].strip() if len(parts) > 1 else name)[:BC_LEN]
            listed[key.encode()] = name.encode()
    if not listed:
        raise ValueError(f'No barcodes in {path}')
    return listed


def build_lookup(listed):
    """Barcode read -> name: exact matches, plus unambiguous 1-mismatch ones."""
    lookup = dict(listed)
    seen = {}
    for key in listed:
        for nb in neighbours(key):
            if nb in listed:
                continue
            seen[nb] = None if nb in seen else key   # ambiguous neighbours are discarded
    lookup.update({nb: listed[key] for nb, key in seen.items() if key is not None})
    return lookup


def window(seq, rule):
    mode, offset = rule
    w = seq[offset:offset + BC_LEN]
    return revcomp(w) if mode == 'rc' else w


def records(handle):
    while True:
        head = handle.readline()
        if not head:
            return
        yield head, handle.readline(), handle.readline(), handle.readline()


def clean_name(head):
    """b'@name/1 comment' -> b'@name' (the mates of a pair must share the name)."""
    name = head.split()[0]
    return name[:-2] if len(name) > 3 and name[-2:-1] == b'/' and name[-1:].isdigit() else name


def choose_rule(listed):
    """Barcode window (orientation, offset) that best matches the cell list in the first reads."""
    if args.rule != 'auto':
        mode, offset = args.rule.split(':')
        return mode, int(offset)

    seqs = []
    with gzip.open(args.r2, 'rb') as fh:
        for rec in itertools.islice(records(fh), args.detect_reads):
            seqs.append(rec[1].strip())
    if not seqs:
        raise ValueError(f'{args.r2} is empty')
    length = len(seqs[0])
    if length < BC_LEN:
        raise ValueError(f'Barcode read of {length} nt, shorter than {BC_LEN}')

    rates = {}
    for mode in ['fwd', 'rc']:
        for offset in range(length - BC_LEN + 1):
            hits = sum(window(s, (mode, offset)) in listed for s in seqs)
            rates[(mode, offset)] = hits / len(seqs)
    best = max(rates, key=rates.get)
    if rates[best] < args.min_match_rate:
        top = sorted(rates.items(), key=lambda x: -x[1])[:3]
        raise ValueError(f'No barcode window matches the cell list (best: {top}). Check the list and its barcode space')
    print(f'Barcode window: {best[0]}:{best[1]} ({rates[best]:.1%} of {len(seqs)} reads match a listed cell)', file=sys.stderr)
    return best


def gunzip(path):
    cmd = [PIGZ, '-dc', path] if PIGZ else ['gzip', '-dc', path]
    return subprocess.Popen(cmd, stdout=subprocess.PIPE, bufsize=BUFSIZE)


def gzip_to(path):
    cmd = [PIGZ, '-3', '-p', str(args.compress_threads), '-c'] if PIGZ else ['gzip', '-3', '-c']
    out = open(path, 'wb')
    return subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=out, bufsize=BUFSIZE), out


def finish(proc, what, out=None):
    """Wait for a gzip process and fail if it did (e.g. truncated input)."""
    if proc.stdin is not None:
        proc.stdin.close()
    if proc.stdout is not None:
        proc.stdout.close()
    if proc.wait() != 0:
        raise RuntimeError(f'gzip failed on {what}')
    if out is not None:
        out.close()


##


def main():

    listed = read_listed(args.cell_barcodes)
    rule = choose_rule(listed)
    lookup = build_lookup(listed)
    mode, offset = rule
    stop = offset + BC_LEN
    # b'_BC:<name>\n' appended to the read name, per barcode
    tags = {key: b'_BC:' + name + b'\n' for key, name in lookup.items()}
    n = kept = exact = 0

    p1, p2, p3 = gunzip(args.r1), gunzip(args.r2), gunzip(args.r3)
    q1, out1 = gzip_to(f'{args.out_prefix}_R1.fastq.gz')
    q3, out3 = gzip_to(f'{args.out_prefix}_R3.fastq.gz')
    f1, f2, f3 = p1.stdout, p2.stdout, p3.stdout
    w1, w3 = q1.stdin, q3.stdin
    get = tags.get
    plus = b'+\n'

    for h1 in f1:
        h2 = f2.readline()
        h3 = f3.readline()
        if not h2 or not h3:
            raise ValueError('R1, R2 and R3 have different numbers of reads')
        s1 = f1.readline()
        f1.readline()
        qual1 = f1.readline()
        seq2 = f2.readline()
        f2.readline()
        f2.readline()
        s3 = f3.readline()
        f3.readline()
        qual3 = f3.readline()
        n += 1
        key = seq2[offset:stop] if mode == 'fwd' else seq2[offset:stop][::-1].translate(COMP)
        tag = get(key)
        if tag is None:
            continue
        kept += 1
        exact += key in listed
        name = clean_name(h1) + tag
        w1.write(name + s1 + plus + qual1)
        w3.write(name + s3 + plus + qual3)

    if f2.readline() or f3.readline():
        raise ValueError('R1, R2 and R3 have different numbers of reads')

    for proc, what in [(p1, args.r1), (p2, args.r2), (p3, args.r3)]:
        finish(proc, what)
    finish(q1, f'{args.out_prefix}_R1.fastq.gz', out1)
    finish(q3, f'{args.out_prefix}_R3.fastq.gz', out3)

    print(f'{n} read pairs, {kept} with a listed cell barcode ({exact} exact)')


##


if __name__ == '__main__':
    main()
