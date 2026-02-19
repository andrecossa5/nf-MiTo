#!/usr/bin/python

# prep_MAESTER script

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='distances',
    description=
    """
    Prepare input for tree building: character/distance matrices and sequences (.fasta) file.
    """
)

# Add arguments
my_parser.add_argument(
    '--afm',
    type=str,
    default='..',
    help='Path to preprocessed afm.h5ad file. Default: .. .'
)

my_parser.add_argument(
    '--n_cores',
    type=int,
    default='8',
    help='n cores to use. Default: 8.'
)

my_parser.add_argument(
    '--boot_strategy',
    type=str,
    default='feature_resampling',
    help='Bootstrap method. Default: feature_resampling.'
)

my_parser.add_argument(
    '--boot_replicate',
    type=str,
    default='original',
    help='Name of the boot replicate. Default: original.'
)

my_parser.add_argument(
    '--frac_char_resampling',
    type=float,
    default=.8,
    help='Fraction of characters to resample at each bootstrapping iteration. Default: .8.'
)


##


# Parse arguments
args = my_parser.parse_args()

path_afm = args.afm
n_cores = args.n_cores
boot_strategy = args.boot_strategy
boot_replicate = args.boot_replicate
frac_char_resampling = args.frac_char_resampling


##


########################################################################

# Preparing run: import code, prepare directories

# Code
import scanpy as sc
import mito as mt
import anndata
anndata.settings.allow_write_nullable_strings = True

########################################################################

# Main
def main():

    # Load afm and extract necessary .uns slots
    afm = sc.read(path_afm)
    metric = afm.uns['distance_calculations']['distances']['metric']
    bin_method = afm.uns['genotyping']['bin_method']
    binarization_kwargs = afm.uns['genotyping']['binarization_kwargs']
    scLT_system = afm.uns['scLT_system']

    # Prep bootstrap kwargs
    kwargs = {
        'boot_replicate':boot_replicate,
        'boot_strategy':boot_strategy,
        'frac_char_resampling':frac_char_resampling
    }

    # Bootstrap
    if scLT_system in ['MAESTER', 'RedeeM']:
        afm_new = mt.tl.bootstrap_MiTo(afm, **kwargs)
    elif scLT_system in ['scWGS', 'Cas9', 'EPI-clone']:
        afm_new = mt.tl.bootstrap_bin(afm, **kwargs)
    else:
        raise ValueError(f'{scLT_system} unavailable! Choose a scLT_system among MAESTER, RedeeeM, scWGS and Cas9, EPI-clone.')

    # Compute distances, recalling
    mt.pp.compute_distances(
        afm_new,
        metric=metric,
        precomputed=False if scLT_system in ['MAESTER', 'RedeeM'] else True, # Bin layer after bootstrapping
        bin_method=bin_method,
        binarization_kwargs=binarization_kwargs,
        ncores=n_cores
    )

    # Save bootstrapped matrix
    afm_new.write('afm_new.h5ad')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
