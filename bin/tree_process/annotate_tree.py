#!/usr/bin/python

# Build cassiopeia

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='annotate_tree',
    description="Process and annotate tree."
)

# Add arguments
my_parser.add_argument(
    '--afm',
    type=str,
    default=None,
    help='Path to afm in afm.h5ad format. Default: .. .'
)

my_parser.add_argument(
    '--tree',
    type=str,
    default=None,
    help='Path to tree in .newick format. Default: .. .'
)

my_parser.add_argument(
    '--annotate_tree',
    type=str,
    default=None,
    help='Path to tree in .newick format. Default: None.'
)

my_parser.add_argument(
    '--max_fraction_unassigned',
    type=float,
    default=.05,
    help='Max fraction of unassigned cells. Default: .05.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_tree = args.tree
path_afm = args.afm
annotate_tree = args.annotate_tree

##


########################################################################

# Preparing run: import code, prepare directories, set logger

# Code
import pickle
import scanpy as sc
import pandas as pd
import mito as mt
import anndata
anndata.settings.allow_write_nullable_strings = True


##


########################################################################

# Main
def main():

    # Prep annot
    afm = sc.read(path_afm)
    cell_meta = afm.obs
    X_raw = pd.DataFrame(afm.X.toarray(), index=afm.obs_names, columns=afm.var_names)
    X_bin = pd.DataFrame(afm.layers['bin'].toarray(), index=afm.obs_names, columns=afm.var_names)
    D = pd.DataFrame(afm.obsp['distances'].toarray(), index=afm.obs_names, columns=afm.obs_names)

    # Load tree
    tree = mt.io.read_newick(path_tree, X_raw=X_raw, X_bin=X_bin, D=D, meta=cell_meta)

    # Cut and annotate tree
    if annotate_tree == "MiTo":
        model = mt.tl.MiToTreeAnnotator(tree)
        model.clonal_inference(max_fraction_unassigned=args.max_fraction_unassigned)

    # Write as pickle
    with open('annotated_tree.pickle', 'wb') as f:
        pickle.dump(tree.copy(), f)


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
