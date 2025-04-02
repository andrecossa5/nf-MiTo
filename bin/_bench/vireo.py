#!/usr/bin/python

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import pickle
from kneed import KneeLocator
from vireoSNP import BinomMixtureVB
import mito as mt
from sklearn.metrics import normalized_mutual_info_score


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='vireo',
    description=
    """
    vireoSNP clonal inference, bayesian (ELBO-based) hyper-parameter 
    tuning.
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
    d['method'] = 'vireoSNP'

    # Read matrix
    afm = sc.read(args.path_afm)

    # Here we go
    _ELBO_mat = []
    maxK = args.maxK
    for k in range(2,maxK+1):
        print(f'Clone n: {k}')
        _model = BinomMixtureVB(n_var=afm.shape[1], n_cell=afm.shape[0], n_donor=k)
        _model.fit(afm.layers['AD'].T, afm.layers['site_coverage'].T, min_iter=30, max_iter=500, max_iter_pre=250, n_init=50, random_seed=1234)
        _ELBO_mat.append(_model.ELBO_inits)

    _ELBO_mat = np.row_stack(_ELBO_mat)
    x = range(2,maxK+1)
    y = np.median(_ELBO_mat, axis=1)
    knee = KneeLocator(x, y).find_knee()[0]
    n_clones = knee

    _model = BinomMixtureVB(n_var=afm.shape[1], n_cell=afm.shape[0], n_donor=n_clones)
    _model.fit(afm.layers['AD'].T, afm.layers['site_coverage'].T, 
               min_iter=30, n_init=50, max_iter=500, max_iter_pre=250, random_seed=1234)

    clonal_assignment = _model.ID_prob
    idx = clonal_assignment.argmax(axis=1)
    labels = np.zeros(idx.shape)
    for i,clone_idx in enumerate(idx):
        labels[i] = clone_idx if clonal_assignment[i,clone_idx]>=.7 else np.nan

    test = np.isnan(labels)
    d['% unassigned'] = test.sum() / labels.size
    d['ARI'] = mt.ut.custom_ARI(afm.obs['GBC'][~test], labels[~test])
    d['NMI'] = normalized_mutual_info_score(afm.obs['GBC'][~test], labels[~test])
    d['labels'] = pd.Series(labels, index=afm.obs_names)   
    d['labels'].loc[lambda x: ~x.isna()] = (
        d['labels']
        .loc[lambda x: ~x.isna()]
        .map(lambda x: f'vireoSNP_{int(x)}')
    )

    # Write
    with open(f'{args.job_id}_vireoSNP.pickle', 'wb') as f:
        pickle.dump(d, f)


# Run
if __name__ == '__main__':
    main()