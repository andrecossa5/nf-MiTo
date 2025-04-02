#!/usr/bin/python

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
import mito as mt
from CCLONE.cluster.NMF import get_wNMF_matrices, NMF_weighted, orth_score
from sklearn.metrics import normalized_mutual_info_score


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='CClone',
    description=
    """
    CClone clonal inference. NMF-based.
    """
)

# Add arguments
my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Sample name. Default: None.'
)
my_parser.add_argument(
    '--job_id', 
    type=str,
    default=None,
    help='job_id. Default: None.'
)
my_parser.add_argument(
    '--path_afm', 
    type=str,
    default='.',
    help='Path to afm.h5ad file. Default: . .'
)
my_parser.add_argument(
    '--maxK', 
    type=int,
    default=15,
    help='Max number of k cluster to try. Default: 15.'
)


# Parse arguments
args = my_parser.parse_args()


##


def main():

    # Dict outputs
    d = {}
    d['sample'] = args.sample
    d['job_id'] = args.job_id
    d['method'] = 'CClone'

    # Read matrix
    afm = sc.read(args.path_afm)

    # Here we go
    afm.layers['ALT'] = afm.layers['AD'].A
    afm.layers['REF'] = afm.layers['site_coverage'].A - afm.layers['AD'].A
    afm.X = afm.X.A

    maxK = args.maxK
    scores = []
    for k in range(2,maxK+1):

        print(f'k={k}')
        afm = get_wNMF_matrices(afm, mode='10X')
        C, _ = NMF_weighted(afm, k=k, max_cycles=100, parallel=True, n_jobs=-1)
        labels = np.argmax(C, axis=1)
        scores.append(orth_score(C))

    k = list(range(2,maxK+1))[np.argmin(scores)]
    C, _ = NMF_weighted(afm, k=k, max_cycles=100, parallel=True, n_jobs=-1)
    labels = np.argmax(C, axis=1)
    test = np.isnan(labels)
    d['% unassigned'] = test.sum() / labels.size
    d['ARI'] = mt.ut.custom_ARI(afm.obs['GBC'][~test], labels[~test])
    d['NMI'] = normalized_mutual_info_score(afm.obs['GBC'][~test], labels[~test])  
    d['labels'] = pd.Series([ f'CClone_{x}' for x in labels ], index=afm.obs_names)

    # Write
    with open(f'{args.job_id}_CClone.pickle', 'wb') as f:
        pickle.dump(d, f)


##


# Run
if __name__ == '__main__':
    main()