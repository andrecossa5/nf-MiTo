# INFER, other scLT tutorial

Despite being optimized fot MT-based scLT (and the MAESTER protocol in particular), nf-MiTo can be used to infer lineages from a variety of other (pre-processed) character matrices.

In this tutorial, we will reproduce the steps needed to infer lineages from character matrices of 3 different scLT systems:

1. The [KP-tracer](10.1016/j.cell.2022.04.015) mouse data, a Cas9-based lineage recorder system
2. The single-cell Whole Genome Sequencing (scWGS) data from [Fabre et al.](https://doi.org/10.1038/s41586-022-04785-z)
3. The [RedeeM](https://doi.org/10.1038/s41586-024-07066-z) data (MT-SNVs enriched from 10x Multiome scATAC-seq libraries)
4. The [EPI-clone](https://doi.org/10.1038/s41586-025-09041-8) data (Tapestri-enriched DNA-methylation status at single-CpG resolution)

**N.B** Here, starting from publicly available resources, we will download, pre-process and build Allele Frequency Matrices (AFMs) ready for nf-MiTo analysis. This analysis is the same described in our [pre-print](https://doi.org/10.1101/2025.06.17.660165) (Figure 5). However, for practical considerations, here we will be using only 1 representative sample for scLT system. The user is encouraged to go through the original works for information about specific scLT systems and pre-processing operations needed. 

# Prerequisites
TO build AFMs, we will need utils from [MiTo](https://github.com/andrecossa5/MiTo). See MiTo's documentation for packge installation and environment setup.

# Cas9-data
To analyze Cas9-based scLT data, we begin collecting KP-Tracer data from [here](https://www.sc-best-practices.org/trajectories/lineage_tracing.html).

```bash
wget "https://zenodo.org/record/5847462/files/KPTracer-Data.tar.gz?download=1"
tar -xvzf KPTracer-Data.tar.gz?download=1
```

Before going further, make sure the folder contains the `KPTracer.alleleTable.FINAL.txt.gz` table.

```bash
ls KPTracer-Data/KPTracer.alleleTable.FINAL.txt.gz
```

From this table, we can now use MiTo I/O functions to generate a AFMs in AnnData format. In this case we will do it for sample 3726_NT_T1.

```python
import mito as mt

# Make AFM
path_ch_matrix = '<your path to KPTracer.alleleTable.FINAL.txt.gz>'
sample = '3726_NT_T1'
afm = mt.io.make_afm(path_ch_matrix, sample=sample, scLT_system='Cas9')
afm

# Write out 
# afm.write('afm.h5ad')
```

That's it. We can now prepare our `--afm_input` as per the [INFER tutorial](2.INFER.md):

<table>
<thead>
<tr>
<th align="left">job_id</th>
<th align="left">sample</th>
<th align="left">afm</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>job_1</code></td>
<td align="left"><code>3726_NT_T1<code></td>
<td align="left"><code>path to 3726_NT_T1 AFM</code></td>
</tr>
</tbody>
</table>

We would also need to prep our `user_params.json` file. In this case:

```json
{  
    "afm_input" : "<path to --afm_input .csv file>",
    "output_folder" : "<path to output folder>",
    "scLT_system" : "Cas9",
    "distance_metric" : "weighted_hamming",
    "annotate_tree" : null
}
```

we set `scLT_system` to `Cas9`, use a suitable distance for this lineage tracer (i.e., `distance_metric` = `weighted_hamming`) and disable tree annotation (only MT-SNVs-based).

Cell phylogeny inference and evaluation can be fired up with:

```bash
nextflow run main.nf -profile docker,local -params-file user_params.json -entry INFER
```

as in the [INFER tutorial](2.INFER.md).


# scWGS data
To analyze scWGS data, we begin collecting data from [here](https://github.com/margaretefabre/Clonal_dynamics). Here, character matrices can be downloaded for individual samples. As an example, we will download filtered mutation calls for sample `PD34493`.

```bash
wget https://raw.githubusercontent.com/margaretefabre/Clonal_dynamics/main/Phylogenies/Files/PD34493/filtered_muts/filtered_muts_PD34493_standard -O filtered_muts_PD34493_standard.RData 
```

This is an .RData file. We need to extract the corresponding chracter matrix in R:

```R

# Load file
load("filtered_muts_PD34493_standard.RData")

# Extract character matrix
bin_matrix <- filtered_muts$COMB_mats.tree.build$Genotype_bin # Binary genotypes

# Write out as csv
output_path <- "<your path here>"
write.csv(bin_matrix, paste0(output_path, "/filtered_muts_PD34493_standard.csv"), row.names=TRUE)
```

We can now build our AFM with MiTo:

```python
import mito as mt

# Make AFM
path_ch_matrix = '<your path to filtered_muts_PD34493_standard.csv>'
sample = 'PD34493'
afm = mt.io.make_afm(path_ch_matrix, sample=sample, scLT_system='scWGS')
afm

# Write out 
# afm.write('afm.h5ad')
```

We can now prepare our `--afm_input` as before:

<table>
<thead>
<tr>
<th align="left">job_id</th>
<th align="left">sample</th>
<th align="left">afm</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left"><code>job_1</code></td>
<td align="left"><code>PD34493<code></td>
<td align="left"><code>path to PD34493 AFM</code></td>
</tr>
</tbody>
</table>

We would also need to prep our `user_params.json` file. In this case:

```json
{  
    "afm_input" : "<path to --afm_input .csv file>",
    "output_folder" : "<path to output folder>",
    "scLT_system" : "scWGS",
    "distance_metric" : "jaccard",
    "annotate_tree" : null
}
```

we set `scLT_system` to `scWGS`, use a suitable distance for this lineage tracer (i.e., `distance_metric` = `jaccard`) and disable tree annotation (only MT-SNVs-based).

Cell phylogeny inference and evaluation can be fired up with:

```bash
nextflow run main.nf -profile docker,local -params-file user_params.json -entry INFER
```

as before.

# RedeeM data
For Redeem data, we will follow the recipe from [Lareau et al.,](https://github.com/caleblareau/redeem-reanalysis-reproducibility), 
for the sample `Young1.T1.HSC.Consensus.final`.

First, we download the data, and inspect its content:

```bash
wget https://figshare.com/ndownloader/files/43428123
mv 43428123 Young1.T1.HSC.Consensus.final.tar.gz
tar -xf Young1.T1.HSC.Consensus.final.tar.gz
tree Young1.T1.HSC.Consensus.final
```

The data should look like this (see [RedeeM-V](https://github.com/chenweng1991/redeemV) for any additional pre-processing information):

```
Young1.T1.HSC.Consensus.final
├── chrM_refAllele.txt
├── QualifiedTotalCts.gz
├── RawGenotypes.Sensitive
├── RawGenotypes.Sensitive.StrandBalance
├── RawGenotypes.Specific
├── RawGenotypes.Specific.StrandBalance
├── RawGenotypes.Total
├── RawGenotypes.Total.StrandBalance
├── RawGenotypes.VerySensitive
├── RawGenotypes.VerySensitive.StrandBalance
├── RawGenotypes.VerySensitive.StrandBalance_small
├── StrandBiaseBlackList
├── StrandBiase.png
├── TotalRawBamRows.Mito
```

From this variant calls, we will perform an interactive variant filtering, as suggested in [Lareau et al.,](https://github.com/caleblareau/redeem-reanalysis-reproducibility) (see also the ongoing debate on RedeeM data quality).

```R
library(tidyverse)
library(redeemR)

path_input <- '<your path here/Young1.T1.HSC.Consensus.final/'

# Create trimmed-ends basecall table
VariantsGTSummary <- redeemR.read.trim(
    path_input, thr="S", Processed=F, rdsname="/VariantsGTSummary.S.trim5_binom.rds", edge_trim=5
)

# Create filtered redeemR object from raw (trimmed) basecall table
redeemR <- Create_redeemR_model(VariantsGTSummary, qualifiedCellCut=10, VAFcut=1, Cellcut=2)

# Filter redeem calls, last utils from Weng et al. to do that
redeemR <- clean_redeem(redeemR, fdr=0.05)
redeemR <- clean_redeem_removehot(redeemR)

# Extract filtered_basecalls, cell_meta an var_meta tables to build the AFM in python.
write.csv(redeemR@GTsummary.filtered, paste0(path_input, 'filtered_basecalls.csv'))
write.csv(redeemR@CellMeta, paste0(path_input, 'cell_meta.csv'))
write.csv(redeemR@V.fitered,  paste0(path_input, 'var_meta.csv'))
```

With these tables, we can now build our AFM with `MiTo`:

```python
import mito as mt

# Make AFM
path_ch_matrix = '<your path to Young1.T1.HSC.Consensus.final folder>'
sample = 'Young1.T1.HSC.Consensus.final'
afm = mt.io.make_afm(path_ch_matrix, sample=sample, scLT_system='RedeeM')
afm

# Write out 
# afm.write('afm.h5ad')
```

nf-MiTo INFER can be launched on this AFM as previously shown for Cas9 and scWGS matrices, this time with `--scLT_system` == `RedeeM` and `--distance_metric` == `jaccard` or `weighted_jaccard`.


# EPI-clone data
Finally, we can explore EPI-clone methylation data. First, we download the data from mouse experiments (see 
[EPI-clone](https://doi.org/10.1038/s41586-025-09041-8)), and inspect its content:

```bash
wget https://figshare.com/ndownloader/files/42479346 
mv seurat_for_figshare.rds seurat_42479346.rds
```

We can extract from the seurat object the needed slots (i.e., cell metadata and binary methylation profiles):

```R
seurat <- readRDS('seurat_42479346.rds')
write.csv(seurat@meta.data, 'cell_meta.csv')   
write.csv(as.data.frame(t(seurat@assays$DNAm@data)), 'DNAm_binary.csv')
```

We will need also epimutation annotations and selection tables. For this experiment (see pubblication main Fig.1):

```bash
wget https://raw.githubusercontent.com/veltenlab/EPI-clone/main/infos/panel_info_dropout_pwm.tsv
wget https://raw.githubusercontent.com/veltenlab/EPI-clone/main/infos/cpg_selection.csv
```

Our input folder should look something like this:

```
├── DNAm_binary.csv
├── cells_meta.csv
├── cpg_selection.csv
├── panel_info_dropout_pwm.tsv
└── seurat_42479346.rds
```

To build EPI-clones AFMs, we proceed as the usual (N.B. in this case we need to extract cells from single samples, 
so we need to provide also `path_meta` to the `mt.io.make_afm` function).

```python
import mito as mt

# Make AFM
path_ch_matrix = '<your path to EPI-clone data folder>'
path_meta = '<your path to EPI-clone data folder>/cells_meta.csv'
sample = 'LARRY_mouse1'
afm = mt.io.make_afm(path_ch_matrix, sample=sample, path_meta=path_meta, scLT_system='EPI-clone')
afm

# Write out 
# afm.write('afm.h5ad')
```

nf-MiTo INFER can be launched as shown before, with `--scLT_system` == `EPI-clone`. The preferred metric for this scLT lineage tracer is `--distance_metric` == `jaccard`.


# TODO:
* Add functionalities to annotate tree with other than MT-SNVs data (in progress)
* Add workflows for scWGS, Cas9, RedeeM, and EPI-clone pre-processing from raw sequencing data (long term) 



