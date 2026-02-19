#!/usr/bin/python

# onesample script

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='onesample',
    description=
    """
    One sample script: fast afm preprocessing, lineage inference and reporting.
    """
)

# Add arguments

my_parser.add_argument(
    '--path_afm',
    type=str,
    default='..',
    help='Path to <name>.h5ad file from mito_preprocessing pipeline. Default: .. .'
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
    help='job id. Default: None.'
)

my_parser.add_argument(
    '--cell_filter',
    type=str,
    default=None,
    help='Cell filtering method. Default: None.'
)

my_parser.add_argument(
    '--mean_cov_all',
    type=float,
    default=10,
    help='Mean MT-genome coverage. Default: 10.'
)

my_parser.add_argument(
    '--median_cov_target',
    type=float,
    default=25,
    help='Median coverage in MAESTER target regions. Default: 25.'
)

my_parser.add_argument(
    '--min_perc_covered_sites',
    type=float,
    default=.75,
    help='Minimum percentage of covered target sites. Default: .75.'
)

my_parser.add_argument(
    '--filtering',
    type=str,
    default=None,
    help='Variant filtering method. Default: None.'
)

my_parser.add_argument(
    '--min_cell_number',
    type=int,
    default=0,
    help='Min number of cell in <lineage_column> categories to retain them. Default: 1.'
)

my_parser.add_argument(
    '--min_var_quality',
    type=int,
    default=30,
    help='Min Q30 phred score of a MT-SNV ADs. Default: 30.'
)

my_parser.add_argument(
    '--min_cov',
    type=int,
    default=5,
    help='Minimum site coverage. Default: 10.'
)

my_parser.add_argument(
    '--min_frac_negative',
    type=float,
    default=.2,
    help='Minimum fraction of -cells to consider a MT-SNV. Default: .2.'
)

my_parser.add_argument(
    '--min_n_positive',
    type=int,
    default=5,
    help='Minimum number of +cells to consider a MT-SNV. Default: 5.'
)

my_parser.add_argument(
    '--af_confident_detection',
    type=float,
    default=.02,
    help='Allelic Frequency of confident detection. Default: .02.'
)

my_parser.add_argument(
    '--min_n_confidently_detected',
    type=int,
    default=2,
    help='Minimum number of confidently detected +cells to consider a MT-SNV. Default: 2.'
)

my_parser.add_argument(
    '--min_mean_AD_in_positives',
    type=float,
    default=1.25,
    help='Minimum number of mean AD in +cells to consider a MT-SNV. Default: 1.5.'
)

my_parser.add_argument(
    '--min_mean_DP_in_positives',
    type=float,
    default=25,
    help='Minimum number of mean DP in +cells to consider a MT-SNV. Default: 20.'
)

my_parser.add_argument(
    '--t_prob',
    type=float,
    default=.7,
    help='Probability threshold for assigning cells to 0/1 mixture binomial components if bin_method=MI_TO. Default: .7.'
)

my_parser.add_argument(
    '--t_vanilla',
    type=float,
    default=0,
    help='AF threshold to assigning cells to 0/1 genotypes if bin_method=MI_TO or vanilla. Default: 0.'
)

my_parser.add_argument(
    '--min_AD',
    type=int,
    default=1,
    help='Min number of AD to assign a 0/1 genotype. Default: 1.'
)

my_parser.add_argument(
    '--bin_method',
    type=str,
    default='MiTo',
    help='Binarization method. Default: MiTo.'
)

my_parser.add_argument(
    '--k',
    type=int,
    default=5,
    help='k neighbors to smooth genotipe if bin_method==MiTo_smooth. Default: 5.'
)

my_parser.add_argument(
    '--gamma',
    type=float,
    default=.25,
    help='% posterior probability that is smoothed in bin_method==MiTo_smooth. Default: .25.'
)

my_parser.add_argument(
    '--min_n_var',
    type=int,
    default=1,
    help='Min n variants. Default: 1.'
)

my_parser.add_argument(
    '--metric',
    type=str,
    default='weighted_jaccard',
    help='Distance metric for cell-cell distance computation. Default: weighted_jaccard.'
)

my_parser.add_argument(
    '--solver',
    type=str,
    default='UPMGA',
    help='Tree-building algorithm for fast lineage inference. Default: UPMGA.'
)

my_parser.add_argument(
    '--min_cell_prevalence',
    type=float,
    default=.1,
    help='Cell prevalence to use binomial modeling in MiTo genotyping. Default: .1.'
)

my_parser.add_argument(
    '--lineage_column',
    type=str,
    default=None,
    help='Lineage column (i.e., GBC for benchmarks). Default: None.'
)

my_parser.add_argument(
    '--ncores',
    type=int,
    default=1,
    help='n cores for computing distances. Default: 1.'
)

my_parser.add_argument(
    '--filter_dbs',
    type=str,
    default="true",
    help='Filter MT-SNVs with dbSNP and REDIdb database. Default: true.'
)

my_parser.add_argument(
    '--filter_moran',
    type=str,
    default="true",
    help='Filter MT-SNVs with spatial auto-correlation. Default: true.'
)

my_parser.add_argument(
    '--max_fraction_unassigned',
    type=float,
    default=.1,
    help='Max fraction of unassigned cells. Default: .1.'
)

my_parser.add_argument(
    '--spatial_metrics',
    type=str,
    default="true",
    help='Add spatial metrics. Default: true.'
)


##


# Parse arguments
args = my_parser.parse_args()


##


########################################################################

# Preparing run: import code, prepare directories

# Code
import os
import numpy as np
import pandas as pd
import scanpy as sc
import mito as mt
import anndata
anndata.settings.allow_write_nullable_strings = True

########################################################################


def main():

    # Extract kwargs
    cell_filter, kwargs, filtering_kwargs, \
    binarization_kwargs, tree_kwargs = mt.ut.extract_kwargs(args)

    # Filter matrix and calculate metrics
    afm = sc.read(args.path_afm)
    afm = mt.pp.filter_cells(
        afm,
        cell_filter=cell_filter,
        mean_cov_all=args.mean_cov_all,
        median_cov_target=args.median_cov_target,
        min_perc_covered_sites=args.min_perc_covered_sites
    )
    afm, tree = mt.pp.filter_afm(
        afm,
        filtering_kwargs=filtering_kwargs,
        binarization_kwargs=binarization_kwargs,
        tree_kwargs=tree_kwargs,
        return_tree=True,
        **kwargs
    )

    # MiTo clones
    model = mt.tl.MiToTreeAnnotator(tree)
    model.clonal_inference(max_fraction_unassigned=args.max_fraction_unassigned)

    # Save info

    # Options
    options = {}
    options['scLT_system'] = afm.uns['scLT_system']
    options['pp_method'] = afm.uns['pp_method']
    options['min_cell_number'] = kwargs['min_cell_number']
    options['lineage_column'] = kwargs['lineage_column']
    options['filtering'] = kwargs['filtering']
    options['bin_method'] = kwargs['bin_method']
    options['min_n_var'] = kwargs['min_n_var']
    options['filter_dbs'] = kwargs['filter_dbs']
    options['filter_moran'] = kwargs['filter_moran']
    options = {
        **options, **filtering_kwargs, **binarization_kwargs,
        **afm.uns['cell_filter'], **tree_kwargs
    }
    (
        pd.Series(options)
        .to_frame('value').reset_index(names='option')
        .assign(sample=args.sample, job_id=args.job_id)
        [['sample', 'job_id', 'option', 'value']]
        .to_csv(f'job_{args.job_id}_options.csv', index=False, header=False)
    )

    # Metrics
    tree = model.tree.copy()
    metrics = {}
    metrics['unassigned'] = tree.cell_meta['MiTo clone'].isna().sum()
    metrics['n MiTo clone'] = tree.cell_meta['MiTo clone'].nunique()
    metrics['corr'] = mt.ut.calculate_corr_distances(tree)[0]
    metrics['mean_CI'] = np.median(mt.ut.CI(tree))
    metrics['mean_RI'] = np.median(mt.ut.RI(tree))

    lineage_column = kwargs['lineage_column']
    if lineage_column is not None and lineage_column in tree.cell_meta.columns:
        metrics[f'n_{lineage_column}_groups'] = tree.cell_meta[lineage_column].nunique()
        metrics['AUPRC'] = mt.ut.distance_AUPRC(afm.obsp['distances'].toarray(), afm.obs[lineage_column])
        test = tree.cell_meta['MiTo clone'].isna()
        metrics['ARI'] = mt.ut.custom_ARI(tree.cell_meta.loc[~test, lineage_column], tree.cell_meta.loc[~test, 'MiTo clone'])
        metrics['NMI'] = mt.ut.normalized_mutual_info_score(
            tree.cell_meta.loc[~test, lineage_column], tree.cell_meta.loc[~test, 'MiTo clone']
        )

    if 'raw_basecalls_metrics' in afm.uns:
        metrics.update(afm.uns['raw_basecalls_metrics'])

    metrics.update(afm.uns['dataset_metrics'])
    metrics['n_dbSNP'] = afm.uns['char_filter']['n_dbSNP']
    metrics['n_REDIdb'] = afm.uns['char_filter']['n_REDIdb']
    metrics['n_not_autocorrelated'] = afm.uns['char_filter']['n_not_autocorrelated']

    # Save
    (
        pd.Series(metrics)
        .to_frame('value').reset_index(names='metric')
        .assign(sample=args.sample, job_id=args.job_id)
        [['sample', 'job_id', 'metric', 'value']]
        .to_csv(f'job_{args.job_id}_metrics.csv', index=False, header=False)
    )


##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
