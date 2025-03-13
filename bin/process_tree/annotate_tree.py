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


##


########################################################################

# Main
def main():

    # Prep annot
    afm = sc.read(path_afm)
    cell_meta = afm.obs
    X_raw = pd.DataFrame(afm.X.A, index=afm.obs_names, columns=afm.var_names)
    X_bin = pd.DataFrame(afm.layers['bin'].A, index=afm.obs_names, columns=afm.var_names)
    D = pd.DataFrame(afm.obsp['distances'].A, index=afm.obs_names, columns=afm.obs_names)

    # # Filter MT-SNVs, as in build_tree
    # if afm.uns['scLT_system'] != 'Cas9':
    #     test_not_germline = ((X_bin==1).sum(axis=0) / X_bin.shape[0]) <= .95        # Prevalence <= 95%
    #     test_not_too_rare = (X_bin==1).sum(axis=0) >= 2                             # At least 2 cells
    #     test = (test_not_germline) & (test_not_too_rare)
    #     X_raw = X_raw.loc[:,test].copy()
    #     X_bin = X_bin.loc[:,test].copy()

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