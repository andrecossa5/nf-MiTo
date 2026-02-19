#!/usr/bin/python

# prep_MAESTER script

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='distance_metrics',
    description=
    """
    Prepare input for tree building: character/distance matrices and sequences (.fasta) file.
    """
)

# Add arguments
my_parser.add_argument(
    '--afm',
    type=str,
    default=None,
    help='List of paths to all bootstrapped afm.h5ad files. Default: None'
)

my_parser.add_argument(
    '--replicates',
    type=str,
    default=None,
    help='List of replicates identifiers mapping to all afm.h5ad files. Default: None'
)

my_parser.add_argument(
    '--job_id',
    type=str,
    default=None,
    help='job_id identifier. Default: None .'
)

my_parser.add_argument(
    '--n_cores',
    type=int,
    default=8,
    help='n cores to use. Default: 8.'
)

my_parser.add_argument(
    '--K',
    type=int,
    default=10,
    help='k NN for NN metrics. Default: 10.'
)

my_parser.add_argument(
    '--lineage_column',
    type=str,
    default=None,
    help='Lineage column for benchmark. Default: None.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_afm = args.afm.strip('[|]').split(', ')
replicates = args.replicates.strip('[|]').split(', ')
K = args.K
job_id = args.job_id
lineage_column = args.lineage_column


##


########################################################################

# Preparing run: import code, prepare directories

# Code
import numpy as np
import scanpy as sc
import mito as mt
import anndata
anndata.settings.allow_write_nullable_strings = True

########################################################################

# Main
def main():

    # Metrics d
    metrics = {}

    # Load distance matrices
    DISTANCES = {}
    for k,v in zip(replicates, path_afm):
        DISTANCES[k] = sc.read(v).obsp['distances']

    # Load cell meta
    idx_observed = [ i for i,x in enumerate(replicates) if x == 'observed' ][0]
    afm = sc.read(path_afm[idx_observed])
    meta = afm.obs

    if lineage_column is not None and lineage_column in meta.columns:

        labels = meta[lineage_column].astype(str)

        # kNN metrics
        D = DISTANCES['observed'].toarray()
        idx = mt.pp.kNN_graph(D=D, k=K, from_distances=True)[0]
        # _, _, acc_rate = mt.ut.kbet(idx, labels, only_score=False)
        median_entropy = mt.ut.NN_entropy(idx, labels)
        median_purity = mt.ut.NN_purity(idx, labels)
        # metrics['kBET_rejection_rate'] = 1-acc_rate
        metrics['median_NN_entropy'] = median_entropy
        metrics['median_NN_purity'] = median_purity

        # AUPRC
        metrics['AUPRC'] = mt.ut.distance_AUPRC(D, labels)

    # Corr
    L = []
    for k in DISTANCES:
        L.append(DISTANCES[k].toarray().flatten())
    del DISTANCES
    metrics['corr'] = np.mean(np.corrcoef(np.array(L)))

    # Save
    afm.uns['distance_metrics'] = metrics
    afm.write('afm_filtered.h5ad')


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################


