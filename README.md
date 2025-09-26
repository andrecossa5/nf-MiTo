[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-blue.svg)](https://www.docker.com/)

# nf-MiTo

Integrated Nextflow pipeline for mitochondrial SNV-based single-cell lineage tracing and multi-omics.

## Overview

nf-MiTo is a flexible Nextflow pipeline for mitochondrial SNV-based single-cell lineage tracing (MT-scLT). 
The pipeline is specifically focused on [MAESTER](10.1038/s41587-022-01210-8) data, where MT-SNVs are profiled from full lenght 10x cDNA. In particular, nf-MiTo implements end-to-end workflows for coupled analysis of Gene Expression (GEX) and MT-SNVs (MT) data raw sequencing data (.fastq or .bam format):

1. Preprocessing (raw sequencing data)
2. Allele Frequency Matrix (AFM) filtering
3. Cell Genotyping
4. Calculation of cell-cell distances
5. Lineage inference

nf-MiTo supports lineage inference (step 2-5, from pre-processed character matrices) of other lineage tracing systems, including RedeeM (MT-DNA, enriched from 10x multiome libraries), Cas9-based lineage recorders, and single-cell Whole Genome Sequencing data (scWGS), representing a novel unifying framework for scLT data analysis.

See our recent pre-print [MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics](https://doi.org/10.1101/2025.06.17.660165) for detailed explanations and benchmarks.

### Key Features

- 📁 **Multiple raw data input formats**: raw FASTQ files mitochondrial BAM files
- 🧬 **Multiple lineage tracing systems**: MAESTER, RedeeM, Cas9-based, and scWGS systems
- ⚡  **Highly optimized methods for MT-scLT data analysis**: [MiTo](https://github.com/andrecossa5/MiTo) python package functionalities, at scale
- 🌳 **Flexible phylogeny reconstruction**: different phylogeny reconstruction and boostrap algorithms
- 🎯 **Estensive parameter optimization**: built-in parameter tuning workflow
- ** Visual Exploration of alternative character spaces and associated trees.
- ☁️ **Scalable execution**: local, HPC, and cloud environment support

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (≥22.04.0)
- [Docker](https://www.docker.com/), [Singularity](https://sylabs.io/singularity/) (preferred), or [Conda](https://conda.io/)

### Basic Usage

nf-MiTo runs on any machine/HPC cluster supporting Docker/Singularity containers:

```bash
nextflow run main.nf -c <user.config> -params-file <params.json> -profile <chosen_profiles> -entry <chosen_entrypoint>
```

With a single command, the user can specify for custom:

- Run configs: `-c <user.config>`, 
- Parameters: `-params-file <params.json>`, 
- Profiles: `-profile <chosen_profiles>`

and opt for the main pipeline workflow (end-to-end), or one of the 4 alternative entrypoints (`PREPROCESS`,
`TUNE`, `EXPLORE`, `INFER`, with the `-entry <chosen_entrypoint>` option). See [configuration](https://www.nextflow.io/docs/latest/config.html) for details.

## Parameters

nf-MiTo implements n=... parameters controlling one or more of the available entrypoints. 
These parameters are grouped in 8 distinct groups. See `nextflow.config` and `nextflow_schema.json` 
for additional type information. See also [MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics](https://doi.org/10.1101/2025.06.17.660165) main text and Supplementary Informations for detailed benchmarks.

### Input/Output parameters

Input/Output Options control nf-MiTo pipeline I/O operations:

<table>
<thead>
<tr>
<th align="left">Parameter</th>
<th align="left">Description</th>
<th align="left">Default</th>
<th align="left">Type</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="white-space: nowrap;"><code>--raw_data_input_type</code></td>
<td align="left">Input data type [<code>fastq</code>, <code>fastq, MAESTER</code>, <code>mitobam</code>]</td>
<td align="left"><code>mitobam</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--raw_data_input</code></td>
<td align="left">Raw data CSV file, preprocessing workflows</td>
<td align="left"><code>–</code></td>
<td align="left">Path</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--afm_input</code></td>
<td align="left">AFM data CSV file, inference workflows</td>
<td align="left"><code>–</code></td>
<td align="left">Path</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--output_folder</code></td>
<td align="left">Output directory path</td>
<td align="left"><code>–</code></td>
<td align="left">Path</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--path_meta</code></td>
<td align="left">Cell metadata file path</td>
<td align="left"><code>null</code></td>
<td align="left">Path</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--lineage_column</code></td>
<td align="left">Cell metadata column with lineage annotation</td>
<td align="left"><code>null</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--path_tuning</code></td>
<td align="left">Tuning results file from TUNE workflow</td>
<td align="left"><code>null</code></td>
<td align="left">Path</td>
</tr>
</tbody>
</table>


For the default end-to-end workflow or the PREPROCESS entypoint (i.e., `-entry PREPROCESS`), the `--raw_data_input` parameter is required. This parameter points to a CSV file sheet storing IDs and paths
to raw sequencing data. 3 alternative inputs are supported here: 

1. `--raw_data_input_type` = `fastq`. The user need to pre-process both MAESTER and TENX sequening data. MAESTER and TENX FASTQs are passed (`--raw_data_input` parameter) with the following a sample sheet (CSV file):

| sample | fastq_folder | library |
|-----------|-------------|---------|
| `sample_name` | `FASTQs folder path` | `MAESTER` |
| `sample_name` | `FASTQs folder path` | `TENX` |

Each sample is linked to a `fastq_folder` for its `MAESTER` and `TENX` (i.e., GEX) library, respectively. 

2. `--raw_data_input_type` = `fastq, MAESTER`. The user is only interested into MAESTER data pre-processing (i.e., GEX data has been analyzed independently). MAESTER FASTQs are passed (`--raw_data_input` parameter) with the following sample sheet (CSV file):

| sample | fastq_folder | cell_barcodes |
|-----------|-------------|---------|
| `sample_name` | `MT FASTQs folder path` | `cell_barcodes path` |

Each sample is linked to a `fastq_folder` for its `MAESTER` library, and a `cell_barcodes.txt` file with cell barcodes of interests (e.g., qualified cell barcodes from independent GEX data analysis). Only reads associated to these cell barcodes will be filtered during pre-processing of `MAESTER` library.

3. `--raw_data_input_type` = `mitobam`. The user already has aligned MAESTER libraries to the genome reference. These aligned reads (i.e., .bam file) are passed (`--raw_data_input` parameter) with the following a sample sheet (CSV file):

| sample | bam | cell_barcodes |
|-----------|-------------|---------|
| `sample_name` | `MT bam path` | `cell_barcodes path` |

Each sample is linked to a `MT bam file` for its `MAESTER` library, and a `cell_barcodes.txt` file with cell barcodes of interests. The bam file should have `CB` and `UB` tags, specifying for cell barcode and UMI, respectively.

All the other entrypoints (i.e., INFER, TUNE, EXPLORE) do *not* implement raw data pre-processing. Thus, it is assumed that properly formatted AFMs have been generated *before* running these workflows. These AFMs are passed (i.e., `--afm_input`) with the following a sample sheet (CSV file):

| job_id | sample | afm |
|-----------|-------------|---------|
| `job_id` | `sample_name` | `AFM path` |

The `job_id` here serve as label links a certain sample (`sample_name`) with a specific AFM (`AFM path`, e.g. from different raw data pre-processing workflows).

The `--output_folder` parameter is always required from the user.

(Optional) the user can provide its custom cell metadata (`--path_meta` parameter, default is `null`). To be valid, these metadata should be formatted as a CSV file indexed with 10x cell barcodes to which their sample name is appended (`_<sample name>`):

| cell | col1 | col2 | 
|-----------|-------------|---------|
| cell1_sampleA | ... | ... | 
| celln_sampleB | ... | ... | 

`--lineage_column` can specify for an orthogonal lineage measurement (i.e., lentiviral barcodes) or inference (e.g. gene-expression derived cell types). If valid, these labels can be used as to compute variant enrichments in these cell groups, or to quantify the relationship with inferred MT-clones. 

### Reference Genome parameters

Reference Genome parameters control nf-MiTo pipeline resources for alignment of sequencing reads. 

<table>
<thead>
<tr>
<th align="left">Parameter</th>
<th align="left">Description</th>
<th align="left">Default</th>
<th align="left">Type</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="white-space: nowrap;"><code>--build_STAR_index</code></td>
<td align="left">Build STAR index from reference files</td>
<td align="left"><code>true</code></td>
<td align="left">Boolean</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--prebuilt_STAR_index</code></td>
<td align="left">Path to prebuilt STAR index directory</td>
<td align="left"><code>null</code></td>
<td align="left">Path</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--reference_genome</code></td>
<td align="left">Reference genome FASTA file URL or path</td>
<td align="left"><a href="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.p13.genome.fa.gz">GENCODE GRCh38.p13</a></td>
<td align="left">Path/URL</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--gtf</code></td>
<td align="left">GTF annotation file URL or path</td>
<td align="left"><a href="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz">GENCODE v32 GTF</a></td>
<td align="left">Path/URL</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--NUMTs_regions</code></td>
<td align="left">Nuclear mitochondrial DNA regions BED file URL or path</td>
<td align="left"><a href="https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed">hg38 NUMTs blacklist</a></td>
<td align="left">Path/URL</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--tenx_whitelist</code></td>
<td align="left">10x Genomics cell barcode whitelist URL or path</td>
<td align="left"><a href="https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz">10x v3 whitelist</a></td>
<td align="left">Path/URL</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--string_MT</code></td>
<td align="left">Mitochondrial chromosome identifier</td>
<td align="left"><code>chrM</code></td>
<td align="left">String</td>
</tr>
</tbody>
</table>

By default, nf-MiTo pulls a reference genome and its annotations from Gencode, and build a STAR index from scratch (i.e., `--build_STAR_index` = `true`). In alternative custom genomes and annotations can be passed (i.e., downloadable URLs, `--reference_genome` and `--gtf` options), or, one can pre-build its own STAR index, and then pass it with the `--prebuilt_STAR_index` option (see the [PREPROCESS](docs/PREPROCESS.md) and [create_STAR_index](docs/create_STAR_index.md) tutorials). 

### scLT System parameters

scLT system parameters specify for the scLT system at hand, and (if `--scLT_system` = `MAESTER`) the chosen
pre-processing tool.

<table>
<thead>
<tr>
<th align="left">Parameter</th>
<th align="left">Description</th>
<th align="left">Default</th>
<th align="left">Type</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="white-space: nowrap;"><code>--scLT_system</code></td>
<td align="left">scLT system [<code>MAESTER</code>, <code>RedeeM</code>, <code>Cas9</code>, <code>scWGS</code>]</td>
<td align="left"><code>MAESTER</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--pp_method</code></td>
<td align="left">Preprocessing method [<code>magatk</code>, <code>mito_preprocessing</code>, <code>cellsnp-lite</code>, <code>samtools</code>, <code>freebayes</code>]</td>
<td align="left"><code>maegatk</code></td>
<td align="left">String</td>
</tr>
</tbody>
</table>

If `--scLT_system` != `MAESTER`, AFMs needs to be generated with MiTo utilities. See ... vignette.

### Sequencing Data Preprocessing parameters

Sequencing Data Preprocessing parameters control how MT-reads are pre-processed (i.e., main workflow or PREPROCESS entrypoint).

<table>
<thead>
<tr>
<th align="left">Parameter</th>
<th align="left">Description</th>
<th align="left">Default</th>
<th align="left">Type</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="white-space: nowrap;"><code>--CBs_chunk_size</code></td>
<td align="left">Cell barcodes chunk size</td>
<td align="left"><code>3000</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--fgbio_UMI_consensus_mode</code></td>
<td align="left">UMI-based read grouping strategy</td>
<td align="left"><code>Identity</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--fgbio_UMI_consensus_edits</code></td>
<td align="left">Max UMI edit distance for read grouping</td>
<td align="left"><code>0</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--fgbio_min_reads_mito</code></td>
<td align="left">Min reads for consensus sequence generation</td>
<td align="left"><code>3</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--fgbio_base_error_rate_mito</code></td>
<td align="left">Max fraction of discordant bases in a read group</td>
<td align="left"><code>0.25</code></td>
<td align="left">Float</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--fgbio_base_quality</code></td>
<td align="left">Min base-calling quality of consensus bases</td>
<td align="left"><code>30</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--fgbio_min_alignment_quality</code></td>
<td align="left">Minimum alignment quality of consensus reads</td>
<td align="left"><code>60</code></td>
<td align="left">Integer</td>
</tr>
</tbody>
</table>

`--CBs_chunk_size` sets the chunk size (i.e., number of cell barcodes) into which the MT-library 
bam file is splitted for parallel pre-processing. All the other parameters control how single-cell
MT-libraries are processed (N.B. when `--pp_method` is `maegatk` or `mito_preprocessing`, see [MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics](https://doi.org/10.1101/2025.06.17.660165) supplementary information for a detailed breakdown of this parameters and their usage in different pre-processing pipelines). In general, [fgbio](https://fulcrumgenomics.github.io/fgbio/) is used for consensus sequence generation. Resulting consunsus reads are then processed to obtain single-cell allelic tables storing consensus UMI counts for all bases (i.e., basecalls) at each MT-genome position. These basecalls are annotated with different statistics (e.g., average base calling quality and consensus score).  

### Cell Filtering parameters

Cell Filtering parameters parameters specify for cell Quality Control (QC) operations.

<table>
<thead>
<tr>
<th align="left">Parameter</th>
<th align="left">Description</th>
<th align="left">Default</th>
<th align="left">Type</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_nUMIs</code></td>
<td align="left">Minimum UMIs per cell</td>
<td align="left"><code>500</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_n_genes</code></td>
<td align="left">Minimum genes per cell</td>
<td align="left"><code>250</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--max_perc_mt</code></td>
<td align="left">Maximum mitochondrial read percentage</td>
<td align="left"><code>0.15</code></td>
<td align="left">Float</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--n_mads</code></td>
<td align="left">n MADs</td>
<td align="left"><code>3</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--cell_filter</code></td>
<td align="left">MT-library filtering strategy [<code>null</code>, <code>filter1</code>, <code>filter2</code>]</td>
<td align="left"><code>filter2</code></td>
<td align="left">String</td>
</tr>
</tbody>
</table>

`--min_nUMIs`, `--min_n_genes`, `--max_perc_mt` and `--n_mads` are actively used only in the PREPROCESS entrypoint (or main end-to-end workflow), *if* `--raw_input_data_type` = `fastq`. These parameters control how cell barcodes are filtered according to standard QC metrics on their GEX library (i.e., n UMIs, genes and % of UMIs mapping to the MT-genome). Considering the number of UMIs and genes, cells are filtered with adaptive thresholds on the upper bound (`--n_mads` * MADs, Median Absolute Deviations).

On the contrary, the `--cell_filter` parameter is used across all workflows, and specify how cells are filtered according to their MT-library properties (i.e., `--scLT_system` in [`MAESTER`, `RedeeM`]). The `filter2` strategy is optimized for MAESTER libraries, while `filter1` can be used for all MT-protocols. 

### Allele Frequency Matrix Preprocessing parameters

Allele Frequency Matrix Preprocessing parameters control how the AFM generated after raw reads pre-processing
is filtered and how cell genotypes are assigned.

<table>
<thead>
<tr>
<th align="left">Parameter</th>
<th align="left">Description</th>
<th align="left">Default</th>
<th align="left">Type</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="white-space: nowrap;"><code>--filtering</code></td>
<td align="left">MT-SNVs filtering method [<code>null</code>, <code>MiTo</code>, <code>MQuad</code>]</td>
<td align="left"><code>MiTo</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--filter_dbs</code></td>
<td align="left">MT-SNVs database filtering</td>
<td align="left"><code>true</code></td>
<td align="left">Boolean</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--spatial_metrics</code></td>
<td align="left">Spatial metrics calculation</td>
<td align="left"><code>false</code></td>
<td align="left">Boolean</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--filter_moran</code></td>
<td align="left">Spatially segregated MT-SNVs filtering</td>
<td align="left"><code>true</code></td>
<td align="left">Boolean</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_cov</code></td>
<td align="left">Min MT-SNV site coverage</td>
<td align="left"><code>5</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_var_quality</code></td>
<td align="left">Min (average) ALT allele base-calling quality</td>
<td align="left"><code>30</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_frac_negative</code></td>
<td align="left">Min fraction of negative cells</td>
<td align="left"><code>0.2</code></td>
<td align="left">Float</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_n_positive</code></td>
<td align="left">Min n of positive cells</td>
<td align="left"><code>5</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--af_confident_detection</code></td>
<td align="left">AF threshold for "confident" detection</td>
<td align="left"><code>0.02</code></td>
<td align="left">Float</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_n_confidently_detected</code></td>
<td align="left">Min n of confidently detected cells</td>
<td align="left"><code>2</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_mean_AD_in_positives</code></td>
<td align="left">Min mean AD in +cells</td>
<td align="left"><code>1.25</code></td>
<td align="left">Float</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_mean_DP_in_positives</code></td>
<td align="left">Min mean DP in +cells</td>
<td align="left"><code>25</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--t_prob</code></td>
<td align="left">Probability threshold for MiTo genotyping</td>
<td align="left"><code>0.7</code></td>
<td align="left">Float</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_AD</code></td>
<td align="left">Minimum allelic depth</td>
<td align="left"><code>2</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_cell_prevalence</code></td>
<td align="left">Min MT-SNV prevalence for MiTo genotyping</td>
<td align="left"><code>0.05</code></td>
<td align="left">Float</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--t_vanilla</code></td>
<td align="left">AF threshold for vanilla genotyping</td>
<td align="left"><code>0</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--bin_method</code></td>
<td align="left">Genotyping method [<code>MiTo</code>, <code>vanilla</code>]</td>
<td align="left"><code>MiTo</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--min_n_var</code></td>
<td align="left">Minimum number of variants per cell</td>
<td align="left"><code>1</code></td>
<td align="left">Integer</td>
</tr>
</tbody>
</table>

`--min_cov`, `--min_var_quality`, `--min_frac_negative`, `--min_n_positive`, `--af_confident_detection`, `--min_n_confidently_detected`, `--min_mean_AD_in_positives`, `--min_mean_DP_in_positives` are active if `--filtering` = `MiTo`. Otherwise all AFM characters are retained (`--filtering` = `null`) or MT-SNVs filtering is performed with the MQuad method (`--filtering` = `MQuad`). 

### Phylogeny Reconstruction parameters

Phylogeny Reconstruction parameters control phylogeny reconstruction from filtered character matrices.

<table>
<thead>
<tr>
<th align="left">Parameter</th>
<th align="left">Description</th>
<th align="left">Default</th>
<th align="left">Type</th>
</tr>
</thead>
<tbody>
<tr>
<td align="left" style="white-space: nowrap;"><code>--distance_metric</code></td>
<td align="left">Distance metric [<code>weighted_jaccard</code>, <code>jaccard</code>, <code>correlation</code>, <code>cosine</code>]</td>
<td align="left"><code>weighted_jaccard</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--tree_algorithm</code></td>
<td align="left">Algorithm [<code>cassiopeia</code>, <code>iqtree</code>, <code>mpboot</code>]</td>
<td align="left"><code>cassiopeia</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--cassiopeia_solver</code></td>
<td align="left">Solver [<code>UPMGA</code>, <code>NJ</code>, <code>spectral</code>, <code>shared_muts</code>, <code>greedy</code>, <code>max_cut</code>]</td>
<td align="left"><code>UPMGA</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--n_boot_replicates</code></td>
<td align="left">Bootstrap replicates</td>
<td align="left"><code>100</code></td>
<td align="left">Integer</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--boot_strategy</code></td>
<td align="left">Bootstrap strategy [<code>feature_resampling</code>, <code>jacknife</code>]</td>
<td align="left"><code>feature_resampling</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--frac_char_resampling</code></td>
<td align="left">% resampled characters</td>
<td align="left"><code>0.8</code></td>
<td align="left">Float</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--support_method</code></td>
<td align="left">Support method [<code>tbe</code>, <code>fbp</code>]</td>
<td align="left"><code>tbe</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--annotate_tree</code></td>
<td align="left">Tree annotation method</td>
<td align="left"><code>MiTo</code></td>
<td align="left">String</td>
</tr>
<tr>
<td align="left" style="white-space: nowrap;"><code>--max_fraction_unassigned</code></td>
<td align="left">Max fraction of unassigned cells</td>
<td align="left"><code>0.1</code></td>
<td align="left">Float</td>
</tr>
</tbody>
</table>

## Configuration
Any custom configuration can be passed to nf-MiTo with the `-c <user.config>` option. See the [`config/user.config`](config/user.config) file for a minimal example for an HPC environment.

## Examples
Tutorials for the main entrypoints and use-cases can be found [here](docs). 

## Best practices


## Troubleshooting

1. **Memory Errors**: Increase memory allocation in configuration
2. **Input Format Errors**: Validate CSV file structure

### Support

- **Issues**: Report bugs via GitHub Issues
- **Discussions**: Community support via GitHub Discussions

## Citation

If you use nf-MiTo in your research, please cite:

```
MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics. doi:https://doi.org/10.1101/2025.06.17.660165
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
