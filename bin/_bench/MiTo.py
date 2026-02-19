#!/usr/bin/python

import argparse
import os
import scanpy as sc
import pickle
import mito as mt
from sklearn.metrics import normalized_mutual_info_score
import anndata
anndata.settings.allow_write_nullable_strings = True


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='MiTo',
    description=
    """
    MiTo cloanl inference, build-in grid-search.
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


# Parse arguments
args = my_parser.parse_args()


##


def main():

    # Dict outputs
    d = {}
    d['sample'] = args.sample
    d['job_id'] = args.job_id
    d['method'] = 'MiTo'

    # Read matrix
    afm = sc.read(args.path_afm)

    # Here we go
    tree = mt.tl.build_tree(afm, precomputed=True)
    model = mt.tl.MiToTreeAnnotator(tree)
    model.clonal_inference()
    tree = model.tree.copy()
    assert (tree.cell_meta.index == afm.obs_names).all()
    labels = tree.cell_meta['MiTo clone']
    test = tree.cell_meta['MiTo clone'].isna()
    d['labels'] = labels
    d['% unassigned'] = test.sum() / labels.size
    d['ARI'] = mt.ut.custom_ARI(afm.obs['GBC'][~test], labels[~test])
    d['NMI'] = normalized_mutual_info_score(afm.obs['GBC'][~test], labels[~test])

    # Write
    with open(f'{args.job_id}_MiTo.pickle', 'wb') as f:
        pickle.dump(d, f)


# Run
if __name__ == '__main__':
    main()
