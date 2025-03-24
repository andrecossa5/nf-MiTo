#!/usr/bin/python

import sys
import scanpy as sc
import mito as mt 


##


def main():

    path_afm = sys.argv[1]
    afm = sc.read(path_afm)
    d = mt.tl.AFM_to_seqs(afm)
    with open('genotypes.fa', 'w') as f:
        for cell in d:
            f.write(f'>{cell}\n{d[cell]}\n')


if __name__ == '__main__':
    main()
