#!/usr/bin/python

# Build cassiopeia

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='build_cassiopeia',
    description="Build cassiopeia trees"
)

# Add arguments
my_parser.add_argument(
    '--path_afm',
    type=str,
    default='..',
    help='Path to afm.h5ad file. Default: .. .'
)

my_parser.add_argument(
    '--path_tuning',
    type=str,
    default=None,
    help='Path to tuning main folder. Default: None.'
)

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
    help='Job id. Default: None.'
)

my_parser.add_argument(
    '--solver',
    type=str,
    default='UPMGA',
    help='Cassiopeia solver. Default: UPMGA.'
)

my_parser.add_argument(
    '--metric',
    type=str,
    default='weighted_jaccard',
    help='Distance metric. Default: weighted_jaccard.'
)

my_parser.add_argument(
    '--ncores',
    type=int,
    default='8',
    help='n cores to use. Default: 8.'
)

my_parser.add_argument(
    '--boot_replicate',
    type=str,
    default='original',
    help='Name of the boot replicate. Default: original.'
)


##


# Parse arguments
args = my_parser.parse_args()


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

    # Extract kwargs
    _, _, _, _, tree_kwargs = mt.ut.extract_kwargs(args, only_tree=True)

    # Build tree
    afm = sc.read(args.path_afm)
    tree = mt.tl.build_tree(afm, precomputed=True, ncores=args.ncores, **tree_kwargs)

    # Write out
    mt.io.write_newick(tree, path=f'rep_{args.boot_replicate}.newick')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
