#!/usr/bin/python

"""
Collapse the per-chunk mgatk tables of a library into the tables MiTo reads.
Usage: mgatk_collapse.py <mt.fa> <out_dir>   (chunk tables <chunk>.{A,C,G,T,coverage,depth}.txt in cwd)
Writes <out_dir>/{A,C,G,T,coverage}.txt.gz, depth.txt.gz and refAllele.txt
(1-based position <TAB> reference base, as mgatk does).
"""

import sys
import gzip
import shutil
from pathlib import Path


##


fasta = sys.argv[1]
out_dir = Path(sys.argv[2])


##


def main():

    out_dir.mkdir(parents=True, exist_ok=True)

    for name in ['A', 'C', 'G', 'T', 'coverage', 'depth']:
        chunks = sorted(Path('.').glob(f'*.{name}.txt'))
        if not chunks:
            raise FileNotFoundError(f'No chunk tables *.{name}.txt in {Path.cwd()}')
        with gzip.open(out_dir / f'{name}.txt.gz', 'wb', compresslevel=4) as out:
            for chunk in chunks:
                with open(chunk, 'rb') as fh:
                    shutil.copyfileobj(fh, out)

    with open(fasta) as fh:
        seq = ''.join(x.strip() for x in fh if not x.startswith('>')).upper()
    with open(out_dir / 'refAllele.txt', 'w') as fh:
        for pos, base in enumerate(seq, start=1):
            fh.write(f'{pos}\t{base}\n')


##


if __name__ == '__main__':
    main()
