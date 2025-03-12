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

my_parser.add_argument(
    '--tree', 
    type=str,
    default=None,
    help='Path to annotated_tree in .pickle format. Default: .. .'
)

my_parser.add_argument(
    '--lineage_column', 
    type=str,
    default=None,
    help='lineage_column. Default: None .'
)

my_parser.add_argument(
    '--job_id', 
    type=str,
    default=None,
    help='job_id. Default: None.'
)


##


# Parse arguments
args = my_parser.parse_args()
path_tree = args.tree
lineage_column = args.lineage_column if args.lineage_column != 'null' else None
job_id = args.job_id


##


########################################################################

# Code
import pickle
import numpy as np
import pandas as pd
import mito as mt
from itertools import chain

########################################################################

# Main
def main():

    # Metrics
    metrics = {}

    # Load tree
    with open(path_tree, 'rb') as f:
        tree = pickle.load(f)
    
    # n cells
    metrics['n_cells'] = tree.layers['transformed'].shape[0]

    # Character metrics
    n_unique_char_states = np.unique(tree.layers['transformed'].values.flatten()).size
    if n_unique_char_states <=3:      
        # SNVs: RedeeM and MAESTER (binary char matrices: 0,1 states), scWGS (ternary char matrices: 0,0.5,1 states)                                                                                   
        metrics['n_characters'] = tree.layers['transformed'].shape[1]
        metrics['median_char_per_cell'] = np.median(np.sum(tree.layers['transformed']>0, axis=1))
        metrics['density'] = np.sum(tree.layers['transformed'].values>0) / np.product(tree.layers['transformed'].shape)
    else:                                                                   
        # INDELs (-1 state: missing info, 0 INDEL state: uncut, others: unique indel states)
        # see https://www.sc-best-practices.org/trajectories/lineage_tracing.html#
        metrics['n_characters'] = tree.layers['transformed'].apply(lambda x: x[x>0].unique().size, axis=0).sum()
        metrics['median_char_per_cell'] = np.median(np.sum(tree.layers['transformed']>0, axis=1))
        X = tree.layers['transformed'].values.flatten()
        metrics['density'] = 1.0 - (np.sum(X==0) / np.sum(X!=-1))

    # Tree-structure/characters consistency (CI/RI indeces)
    metrics['median_CI'] = np.median(mt.ut.CI(tree))
    metrics['median_RI'] = np.median(mt.ut.RI(tree))

    # Bootstrap support
    df_internal = mt.tl.get_internal_node_stats(tree)
    metrics['median_support'] = df_internal['support'].median()
    if 'MiTo clone' in tree.cell_meta.columns:
        metrics['median_support_mut_clades'] = df_internal.loc[df_internal['mut_clade'],'support'].median()
    max_clade = np.percentile(df_internal['clade_size'], 90)
    metrics['median_support_biggest_clades'] = df_internal.query(f'clade_size>={max_clade}')['support'].median()

    # Time
    metrics['median_time'] = df_internal['time'].median()

    # Corr genetic-tree distances
    corr_dists, p_corr_dists = mt.ut.calculate_corr_distances(tree)
    metrics['corr_distances'] = corr_dists
    metrics['corr_distances_pvalue'] = p_corr_dists

    # Clonal assignment metrics
    if 'MiTo clone' in tree.cell_meta.columns:
        metrics['n_clones'] = tree.cell_meta['MiTo clone'].nunique()
        metrics['median_n_cells_clone'] = tree.cell_meta['MiTo clone'].value_counts().median()
        metrics['min_n_cells_clone'] = tree.cell_meta['MiTo clone'].value_counts().min()
        metrics['max_n_cells_clone'] = tree.cell_meta['MiTo clone'].value_counts().max()
        metrics['median_n_assigned_char_per_clone'] = np.median(
            [ len(x.split(';')) for x in tree.cell_meta['muts'].loc[lambda x: ~x.isna()].unique() ]
        )
        metrics['total_assigned_char'] = len(list(
            chain.from_iterable([ x.split(';') for x in tree.cell_meta['muts'].loc[lambda x: ~x.isna()].unique() ])
        ))

    # Benchmark specifics
    if lineage_column is not None and lineage_column in tree.cell_meta.columns:
        test = ~tree.cell_meta['MiTo clone'].isna()
        metrics['ARI'] = mt.ut.custom_ARI(
            tree.cell_meta.loc[test, lineage_column], tree.cell_meta.loc[test, 'MiTo clone']
        )
        metrics['NMI'] = mt.ut.normalized_mutual_info_score(
            tree.cell_meta.loc[test, lineage_column], tree.cell_meta.loc[test, 'MiTo clone']
        )

    # Collapse tree
    n = len(tree.internal_nodes)
    tree.collapse_mutationless_edges(True)
    n_ = len(tree.internal_nodes)
    metrics['frac_remaining_nodes_after_mutationless_edges_collapse'] = n_/n

    # Save
    to_frame = lambda x: pd.Series(x).to_frame('value').reset_index(names='metric').assign(job_id=job_id)
    to_frame(metrics).to_csv('tree_metrics.csv')


##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################