#!/usr/bin/python

# MAESTER script

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='scWGS',
    description=
    """
    Prepare input for tree building, from (pre-processed) single-colony WGS data Allele Frequency Matrix.
    """
)

# Add arguments

my_parser.add_argument(
    '--path_afm',
    type=str,
    default='.',
    help='Path to afm.h5ad file. Default: . .'
)

my_parser.add_argument(
    '--metric',
    type=str,
    default='jaccard',
    help='Distance metric. Default: jaccard.'
)

my_parser.add_argument(
    '--n_cores',
    type=int,
    default=1,
    help='n cores to use. Default: 1.'
)


##


########################################################################

# Code
import scanpy as sc
import mito as mt
import anndata
anndata.settings.allow_write_nullable_strings = True

########################################################################

# Main
def main():

    # Parse arguments
    args = my_parser.parse_args()

    path_afm = args.path_afm
    metric = args.metric
    n_cores = args.n_cores

    # Calculate distances and write out afm
    afm = sc.read(path_afm)
    afm.uns['genotyping'] = {'layer':'bin', 'bin_method':'', 'binarization_kwargs':{}}
    mt.pp.compute_distances(afm, precomputed=True, metric=metric, ncores=n_cores)
    afm.write('afm.h5ad')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
