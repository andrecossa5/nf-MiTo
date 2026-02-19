#!/usr/bin/python

import argparse
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import pickle
import mito as mt
from sklearn.metrics import silhouette_score, normalized_mutual_info_score
import anndata
anndata.settings.allow_write_nullable_strings = True


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='leiden',
    description=
    """
    Leiden clustering, optimized resolution.
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


##


# Parse arguments
args = my_parser.parse_args()


##


def main():

    # Dict outputs
    d = {}
    d['sample'] = args.sample
    d['job_id'] = args.job_id
    d['method'] = 'leiden'

    # Read matrix
    afm = sc.read(args.path_afm)

    # Here we go
    scores = []
    resolutions = np.linspace(0.5,2.5,50)
    for res in resolutions:
        _, _, conn = mt.pp.kNN_graph(D=afm.obsp['distances'].A, k=15, from_distances=True)
        labels = mt.tl.leiden_clustering(conn, res=res)
        silhouette_avg = silhouette_score(afm.obsp['distances'].A, labels, metric='precomputed')
        scores.append(silhouette_avg)

    labels = mt.tl.leiden_clustering(conn, res=resolutions[np.argmax(scores)])
    d['% unassigned'] = 0
    d['ARI'] = mt.ut.custom_ARI(afm.obs['GBC'], labels)
    d['NMI'] = normalized_mutual_info_score(afm.obs['GBC'], labels)
    d['labels'] = pd.Series([ f'leiden_{x}' for x in labels ], index=afm.obs_names)

    with open(f'{args.job_id}_leiden.pickle', 'wb') as f:
        pickle.dump(d, f)


# Run
if __name__ == '__main__':
    main()
