#!/usr/bin/python

"""
Basic QC.
"""

import os
import argparse


##


# Create the parser
my_parser = argparse.ArgumentParser(
    prog='QC',
    description=
    """
    Basic QC.
    """
)

# Add arguments
my_parser.add_argument(
    '--input', 
    type=str,
    default=None,
    help='Path to filtered/ folder from CellRanger/STARSolo. Default: None'
)

my_parser.add_argument(
    '--path_meta', 
    type=str,
    default=None,
    help='Path to cells_meta.csv. Default: None'
)

my_parser.add_argument(
    '--sample', 
    type=str,
    default=None,
    help='Sample name. Default: None'
)

my_parser.add_argument(
    '--min_nUMIs', 
    type=int,
    default=500,
    help='Min n of UMIs allowed. Default: 500'
)

my_parser.add_argument(
    '--min_n_genes', 
    type=int,
    default=250,
    help='Min n of genes allowed. Default: 250'
)

my_parser.add_argument(
    '--max_perc_mt', 
    type=float,
    default=.15,
    help='Max % of MT-gene counts allowed. Default: .15'
)

my_parser.add_argument(
    '--n_mads', 
    type=int,
    default=5,
    help='n MADs to filter out the higher-end nUMIs/n genes distributions. Default: 5'
)


##


# Parse arguments
args = my_parser.parse_args()
path_input = args.input
path_meta = args.path_meta
sample = args.sample
min_nUMIs = args.min_nUMIs
min_n_genes = args.min_n_genes
max_perc_mt = args.max_perc_mt
n_mads = args.n_mads


##


import os
import numpy as np
import pandas as pd
import scanpy as sc


##


def main():

    # Read data data 
    adata = sc.read_10x_mtx(path_input)
    adata.obs_names = adata.obs_names.map(lambda x: f'{x}_{sample}') # NB: covention needed on meta df

    # Meta or not?
    if path_meta is not None and os.path.exists(path_meta):

        meta = pd.read_csv(path_meta, index_col=0)
        if sample in meta['sample'].unique():
            meta = meta.query('sample==@sample').copy()  
            cell_keep = list(set(meta.index) & set(adata.obs_names))   
            assert len(cell_keep>0)     
            adata.var['n_cells'] = (adata.X>0).sum(axis=0).A1
            adata.var['pct_cells'] = (adata.X>0).sum(axis=0).A1 / adata.shape[0] * 100
            gene_keep = adata.var.loc[lambda x: ~x['mt']|x['ribo']].query('pct_cells>=@min_pct_cell').index
            adata = adata[cell_keep,gene_keep].copy()
            adata.obs = meta.loc[cell_keep]
        else:
            raise ValueError(f'{sample} not in meta! Check path_meta...')

    else:

        # Assign sample name
        adata.obs = adata.obs.assign(sample=sample)
   
        # Basic QC, all samples together (all cells also)
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
        qc_metrics_cells, _ = sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=False)
        adata.obs = adata.obs.join(
            qc_metrics_cells.assign(n_genes=lambda x: (adata.X>0).sum(axis=1).A1)
            [['total_counts', 'n_genes', 'pct_counts_mt', 'pct_counts_ribo']]
        )
        adata.var['n_cells'] = (adata.X>0).sum(axis=0).A1
        adata.var['pct_cells'] = (adata.X>0).sum(axis=0).A1 / adata.shape[0] * 100

        # Filter cells and genes
        mad = lambda x: np.median(np.absolute(x - np.median(x)))
        test_UMIs = (adata.obs['total_counts'] >= min_nUMIs) & \
                    (adata.obs['total_counts'] <= n_mads * mad(adata.obs['total_counts']))
        test_genes = (adata.obs['n_genes'] >= min_n_genes) & \
                     (adata.obs['n_genes'] <= n_mads * mad(adata.obs['n_genes']))
        test_mt = (adata.obs['pct_counts_mt'] <= max_perc_mt)
        test = (test_UMIs) & (test_genes) & (test_mt)
        cell_keep = adata.obs_names[test]
        gene_keep = adata.var.loc[lambda x: ~x['mt']|x['ribo']].query('pct_cells>=.01').index
        adata = adata[cell_keep,gene_keep].copy()

    # Write
    adata.write('adata.h5ad')
    adata.obs_names.to_series().to_csv('cell_barcodes.txt', index=False)


##


# Run
if __name__ == '__main__':
    main()

