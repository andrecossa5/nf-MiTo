#!/usr/bin/python

# MAESTER script

########################################################################

# Libraries
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='MiTo',
    description=
    """
    Prepare input for tree building, from MAESTER or RedeeM Allele Frequency Matrix.
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
    '--path_tuning',
    type=str,
    default=None,
    help='Path to tuning output folder. Default: None.'
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
    help='Variant filtering method. Default: None (i.e., all MT-SNVs variants will be retained).'
)

my_parser.add_argument(
    '--min_cell_number',
    type=int,
    default=0,
    help='Min number of cell in <lineage_column> categories to retain them. Default: 0.'
)

my_parser.add_argument(
    '--min_cov',
    type=int,
    default=5,
    help='Minimum mean coverage of the candidate variant site. Default: 5.'
)

my_parser.add_argument(
    '--min_var_quality',
    type=int,
    default=30,
    help='Min phred score of a MT-SNV ADs. Default: 30.'
)

my_parser.add_argument(
    '--min_frac_negative',
    type=float,
    default=.2,
    help='Minimum fraction of negative (i.e., AF==0) cells to consider a MT-SNV. Default: .2.'
)

my_parser.add_argument(
    '--min_n_positive',
    type=int,
    default=5,
    help='Minimum number of positive (i.e., AF>0) cells to consider a MT-SNV. Default: 5.'
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
    help='Minimum number of confidently detected positive cells to consider a MT-SNV. Default: 2.'
)

my_parser.add_argument(
    '--min_mean_AD_in_positives',
    type=float,
    default=1.25,
    help='Minimum number of mean AD in positive cells to consider a MT-SNV. Default: 1.5.'
)

my_parser.add_argument(
    '--min_mean_DP_in_positives',
    type=float,
    default=25,
    help='Minimum number of mean DP in positive cells to consider a MT-SNV. Default: 20.'
)

my_parser.add_argument(
    '--t_prob',
    type=float,
    default=.7,
    help='Probability threshold for assigning cells to 0/1 mixture binomial components if bin_method=MiTo. Default: .7.'
)

my_parser.add_argument(
    '--t_vanilla',
    type=float,
    default=0,
    help='AF threshold to assigning cells to 0/1 genotypes if bin_method=MiTo or vanilla. Default: 0.'
)

my_parser.add_argument(
    '--min_AD',
    type=int,
    default=1,
    help='Min number of AD to assign a 0/1 genotype. Default: 1.'
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
    help='% posterior probability that is smoothed in bin_method==MiTo_smooth. Default: .2.'
)

my_parser.add_argument(
    '--bin_method',
    type=str,
    default='MiTo',
    help='Binarization method. Default: MiTo.'
)

my_parser.add_argument(
    '--min_n_var',
    type=int,
    default=1,
    help='Min n variants. Default: 1.'
)

my_parser.add_argument(
    '--min_cell_prevalence',
    type=float,
    default=.1,
    help='Min cell prevalence to assign 0/1 genotype with the MiTo method. Default: .1.'
)

my_parser.add_argument(
    '--lineage_column',
    type=str,
    default=None,
    help='Lineage column (i.e., GBC for benchmarks). Default: None.'
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
    default=1,
    help='n cores to use. Default: 1.'
)

my_parser.add_argument(
    '--filter_dbs',
    type=str,
    default="true",
    help='Filter MT-SNVs with dbSNP and REDIdb database. Default: true.'
)

my_parser.add_argument(
    '--spatial_metrics',
    type=str,
    default="false",
    help='Add spatial metrics. Default: false.'
)

my_parser.add_argument(
    '--filter_moran',
    type=str,
    default="true",
    help='Add filtering for spatial autocorrelation. Default: true.'
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
    cell_filter, kwargs, filtering_kwargs, \
    binarization_kwargs, tree_kwargs = mt.ut.extract_kwargs(args)

    # Read AFM
    afm = sc.read(args.path_afm)

    # Filter cells
    if cell_filter is not None:
        afm = mt.pp.filter_cells(
            afm,
            cell_filter=cell_filter,
            mean_cov_all=args.mean_cov_all,
            median_cov_target=args.median_cov_target,
            min_perc_covered_sites=args.min_perc_covered_sites
        )

    # Filter MT-SNVs
    afm = mt.pp.filter_afm(
        afm,
        filtering_kwargs=filtering_kwargs,
        binarization_kwargs=binarization_kwargs,
        tree_kwargs=tree_kwargs,
        return_tree=False,
        **kwargs
    )

    # Write out filtered matrix
    afm.write('afm.h5ad')


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################
