#!/usr/bin/python

# MAESTER script

########################################################################

# Libraries 
import os
import argparse

# Create the parser
my_parser = argparse.ArgumentParser(
    prog='Explore MT-SNVs space',
    description=
    """
    Viz script.
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
    '--covariate', 
    type=str,
    default=None,
    help='Covariate to plot. Default: None.'
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
import pandas as pd
import scanpy as sc
import mito as mt 
import matplotlib.pyplot as plt
import plotting_utils as plu

########################################################################

# Main
def main():

    # Prep folder
    os.mkdir(args.job_id)
    os.chdir(args.job_id)

    # Extract kwargs
    cell_filter, kwargs, filtering_kwargs, \
    binarization_kwargs, tree_kwargs = mt.ut.extract_kwargs(args)

    # Filter cells
    afm_raw = sc.read(args.path_afm)
    afm_raw = mt.pp.filter_cells(afm_raw, cell_filter=cell_filter)
    mt.pp.annotate_vars(afm_raw)
    afm_raw = mt.pp.filter_baseline(afm_raw)

    # Read and format coverage 
    path_coverage = os.path.join(os.path.dirname(args.path_afm), 'tables', 'coverage.txt.gz')
    cov = pd.read_csv(path_coverage, header=None)
    cov.columns = ['pos', 'cell', 'n'] 
    cov['cell'] = cov['cell'].map(lambda x: f'{x}_{args.sample}')
    cov = cov.query('cell in @afm_raw.obs_names')
    cov['cell'] = pd.Categorical(cov['cell'], categories=afm_raw.obs_names)
    cov['pos'] = pd.Categorical(cov['pos'], categories=range(1,16569+1))
    cov = cov.pivot_table(index='cell', columns='pos', values='n', fill_value=0)

    # Filter afm, reduce dimensions
    afm = mt.pp.filter_afm(
        afm_raw,
        filtering_kwargs=filtering_kwargs,
        binarization_kwargs=binarization_kwargs,
        tree_kwargs=tree_kwargs,
        compute_enrichment=True,
        return_tree=False,
       **kwargs
    )
    mt.pp.reduce_dimensions(afm)

    # Build and annotate tree
    tree = mt.tl.build_tree(afm, precomputed=True)
    model = mt.tl.MiToTreeAnnotator(tree)
    model.clonal_inference()


    ##


    # 1. Viz mutation selection
    fig = plt.figure(figsize=(15,4.5))

    ax = fig.add_subplot(1,4,1)
    mt.pl.vars_AF_spectrum(afm_raw, ax=ax, color='#303030', alpha=.7, linewidth=.2)
    mt.pl.vars_AF_spectrum(afm_raw[:,afm.var_names], ax=ax, color='#05A8B3', linewidth=.5, alpha=1)

    ax = fig.add_subplot(1,4,2)
    xticks = [1,2,4,10,30,90,300,1100]
    mt.pl.plot_ncells_nAD(afm_raw, ax=ax,  xticks=xticks, c='#303030', s=2, alpha=.3)
    mt.pl.plot_ncells_nAD(afm, ax=ax, c='#05A8B3', xticks=xticks, s=5, alpha=1, markeredgecolor='k')
    plu.format_ax(ax=ax, ylabel='Mean nAD / +cells', xlabel='n +cells', reduced_spines=True)

    ax = fig.add_subplot(1,4,3, polar=True)
    mt.pl.MT_coverage_polar(cov, var_subset=afm.var_names, ax=ax, 
                      kwargs_subset={'markersize':8, 'c':'#05A8B3'}, 
                      kwargs_main={'c':'#303030', 'linewidth':1.5, 'alpha':.7})

    ax = fig.add_subplot(1,4,4)
    ref_df = mt.ut.load_mt_gene_annot()
    df_plot = ref_df.query('mut in @afm.var_names')
    plu.counts_plot(df_plot, 'Symbol', width=.8, ax=ax, color='#C0C0C0', edgecolor='k', with_label=False)
    plu.format_ax(ax=ax, xticks=df_plot.index, rotx=90, 
                  ylabel='n MT-SNVs', xlabel='Gene', reduced_spines=True)

    fig.subplots_adjust(bottom=.25, top=.8, left=.1, right=.9, wspace=.4)
    fig.savefig('MT_SNVs.png', dpi=500)


    ##


    # 2. Viz mutation profile
    fig = mt.pl.mut_profile(afm.var_names, figsize=(5,2.5))
    fig.tight_layout()
    fig.savefig('mut_profile.png', dpi=500)


    ##


    # 3. Viz embeddings
    cmaps = {
        'MiTo clone' : \
        plu.create_palette(model.tree.cell_meta, 'MiTo clone', sc.pl.palettes.default_102, add_na=True)
    }
    if args.covariate is not None:
        cmaps[args.covariate] = plu.create_palette(
            model.tree.cell_meta, args.covariate, sc.pl.palettes.vega_20_scanpy
        )
        afm.uns[f'{args.covariate}_colors'] = list(cmaps[args.covariate].values())

    fig, ax = plt.subplots(figsize=(9,5))
    sc.pl.embedding(afm, basis='X_umap', color=args.covariate, ax=ax, show=False, save=False, frameon=False)
    fig.subplots_adjust(bottom=.1, top=.9, left=.1, right=.5)
    fig.savefig(f'embeddings.png', dpi=500)


    ##


    # 4. Viz tree
    fig, ax = plt.subplots(figsize=(4.7,5))
    mt.pl.plot_tree(
        tree, ax=ax, 
        features=list(cmaps.keys()), 
        colorstrip_width=5, 
        categorical_cmaps=cmaps,
        feature_internal_nodes='similarity',
        internal_node_subset=model.clonal_nodes,
        internal_node_kwargs={'markersize':8}
    )
    n_clones = tree.cell_meta['MiTo clone'].unique().size
    fig.suptitle(f'n cells: {afm.shape[0]}, n vars: {afm.shape[1]}, n MiTo clones: {n_clones}')
    fig.tight_layout()
    fig.savefig('phylo.png', dpi=500)


    ##


    # 5. Viz tree with muts
    fig, axs = plt.subplots(1,2,figsize=(15,8), gridspec_kw={'wspace': 0.4})

    mt.pl.plot_tree(
        tree, ax=axs[0], 
        colorstrip_spacing=.000001, colorstrip_width=2, 
        orient='down',
        features=list(cmaps.keys()),
        characters=model.ordered_muts, layer='raw',
        categorical_cmaps=cmaps,
        feature_internal_nodes='similarity',
        internal_node_subset=model.clonal_nodes,
        show_internal=True, 
        internal_node_kwargs={'markersize':8}
    )
    plu.add_cbar(
        model.tree.layers['transformed'].values.flatten(), 
        palette='mako', label='AF', 
        ticks_size=8, label_size=9, vmin=.0, vmax=.1,
        ax=axs[0], layout=( (1.02,.3,.02,.2), 'right', 'vertical' )
    )

    mt.pl.plot_tree(
        tree, ax=axs[1],
        colorstrip_spacing=.000001, colorstrip_width=2, 
        orient='down',
        features=list(cmaps.keys()),
        characters=model.ordered_muts, layer='transformed',
        categorical_cmaps=cmaps,
        feature_internal_nodes='similarity',
        internal_node_subset=model.clonal_nodes,
        show_internal=True, 
        internal_node_kwargs={'markersize':8}
    )
    plu.add_legend(
        label='Genotype', ax=axs[1], 
        colors={'REF':'b', 'ALT':'r'}, loc='center left', 
        bbox_to_anchor=(1,.4),
        ticks_size=8, artists_size=10, label_size=9
    )

    fig.tight_layout()
    fig.savefig('phylo_muts.png', dpi=500)


    ##


    # 5. Viz distances
    fig, ax = plt.subplots(figsize=(4.5,4.5))
    mt.pl.heatmap_distances(afm, tree=tree, ax=ax)
    fig.tight_layout()
    fig.savefig('distances.png', dpi=500)


    ##


#######################################################################

# Run program
if __name__ == "__main__":
    main()

#######################################################################

