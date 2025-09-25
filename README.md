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

With a single command, the user can provide its custom:

- Run [configuration](https://www.nextflow.io/docs/latest/config.html): `-c <user.config>`, 
- Parameters [parameters](https://www.nextflow.io/docs/latest/config.html#parameters): `-params-file <params.json>`, 
- Profiles [profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles): `-profile <chosen_profiles>`

and opt for the main pipeline workflow (end-to-end), or one of the 4 alternative entrypoints (`PREPROCESS`,
`TUNE`, `EXPLORE`, `INFER`, with the `-entry <chosen_entrypoint>` option).

## Parameters

nf-MiTo implements n=56 parameters controlling one or more of the available entrypoints.
These parameters are grouped in 8 distinct groups. See also `nextflow.config` and `nextflow_schema.json` 
for additional type information. 

### Input/Output parameters

Input/Output Options control nf-MiTo pipeline I/O operations.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--raw_data_input_type` | Input data type [`fastq`, `fastq, MAESTER`, `mitobam`] | `mitobam` | 
| `--raw_data_input` | Raw data CSV file, preprocessing workflows | `–` |
| `--afm_input` | AFM data CSV file, inference workflows | `–` |
| `--output_folder` | Output directory path | `–` |
| `--path_meta` | Cell metadata file path | `null` |
| `--lineage_column` | Cell metadata column with lineage annotation | `null` |
| `--path_tuning` | Tuning results file from TUNE workflow | `null` |

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

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--build_STAR_index` | Build STAR index from reference files | `true` |
| `--prebuilt_STAR_index` | Path to prebuilt STAR index directory | `null` |
| `--reference_genome` | Reference genome FASTA file URL or path | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh38.p13.genome.fa.gz` |
| `--gtf` | GTF annotation file URL or path | `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.annotation.gtf.gz` |
| `--NUMTs_regions` | Nuclear mitochondrial DNA regions BED file URL or path | `https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed` |
| `--10x_whitelist` | 10x Genomics cell barcode whitelist URL or path | `https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz` |
| `--ref` | Reference genome directory | - |
| `--string_MT` | Mitochondrial chromosome identifier | `chrM` |
| `--whitelist` | 10x v3 whitelist file | - |

The `--ref` and `whitelist` parameters are required only for the main workflow and the PREPROCESS entrypoint.

The `--ref` parameter must pointo to valid STAR index folder, which will be used by [STAR Solo](https://doi.org/10.1101/2021.05.05.442755) for paired-end reads alignment. Following [best practices](10.1038/s41596-022-00795-3) for MT-SNVs data analysis, this STAR index should be build with a reference genome with masked NUMTs. 

The `--whitelist` parameter indicate the 10x v3 cell barcode whitelist, which can be downloaded from [10x](https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist) site. This whitelist is used by STAR Solo to correct cellular barcodes prior than MT-library bam filtering and processing.

### scLT System parameters

scLT system parameters specify for the scLT system at hand, and (if `--scLT_system` = `MAESTER`) the chosen
pre-processing tool.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--scLT_system` | scLT system [`MAESTER`, `RedeeM`, `Cas9`, `scWGS`] | `MAESTER` | |
| `--pp_method` | Preprocessing method [`magatk`, `mito_preprocessing`, `cellsnp-lite`, `samtools`, `freebayes`] | `maegatk` | |

If `--scLT_system` != `MAESTER`, AFMs needs to be generated with MiTo utilities. See [...] vignette.

### Sequencing Data Preprocessing parameters

Sequencing Data Preprocessing parameters control how MT-reads are pre-processed (i.e., main workflow or PREPROCESS entrypoint).

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--CBs_chunk_size` | Cell barcodes chunk size | `3000` |
| `--fgbio_UMI_consensus_mode` | UMI-based read grouping strategy | `Identity` |
| `--fgbio_UMI_consensus_edits` | Max UMI edit distance for read grouping  | `0` |
| `--fgbio_min_reads_mito` | Min reads for consensus sequence generation | `3` |
| `--fgbio_base_error_rate_mito` | Max fraction of discordant bases in a read group | `0.25` |
| `--fgbio_base_quality` | Min base-calling quality of consensus bases | `30` |
| `--fgbio_min_alignment_quality` | Minimum alignment quality of consensus reads | `60` |

`--CBs_chunk_size` sets the chunk size (i.e., number of cell barcodes) into which the MT-library 
bam file is splitted for parallel pre-processing. All the other parameters control how single-cell
MT-libraries are processed (N.B. when `--pp_method` is `maegatk` or `mito_preprocessing`, see [MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics](https://doi.org/10.1101/2025.06.17.660165) supplementary information for a detailed breakdown of this parameters and their usage in different pre-processing pipelines). In general, [fgbio](https://fulcrumgenomics.github.io/fgbio/) is used for consensus sequence generation. Resulting consunsus reads are then processed to obtain single-cell allelic tables storing consensus UMI counts for all bases (i.e., basecalls) at each MT-genome position. These basecalls are annotated with different statistics (e.g., average base calling quality and consensus score).  

### Cell Filtering parameters

Cell Filtering parameters parameters specify for cell Quality Control (QC) operations.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min_nUMIs` | Minimum UMIs per cell | `500` |
| `--min_n_genes` | Minimum genes per cell | `250` |
| `--max_perc_mt` | Maximum mitochondrial read percentage | `0.15` |
| `--n_mads` | n MADs | 3 |
| `--cell_filter` | MT-library filtering strategy [`null`, `filter1`, `filter2`]  | `filter2` |

`--min_nUMIs`, `--min_n_genes`, `--max_perc_mt` and `--n_mads` are actively used only in the PREPROCESS entrypoint (or main end-to-end workflow), *if* `--raw_input_data_type` = `fastq`. These parameters control how cell barcodes are filtered according to standard QC metrics on their GEX library (i.e., n UMIs, genes and % of UMIs mapping to the MT-genome). Considering the number of UMIs and genes, cells are filtered with adaptive thresholds on the upper bound (`--n_mads` * MADs, Median Absolute Deviations).

On the contrary, the `--cell_filter` parameter is used across all workflows, and specify how cells are filtered according to their MT-library properties (i.e., `--scLT_system` in [`MAESTER`, `RedeeM`]). The `filter2` strategy is optimized for MAESTER libraries, while `filter1` can be used for all MT-protocols. See [MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics](https://doi.org/10.1101/2025.06.17.660165) Supplementary Informations for details.

### Allele Frequency Matrix Preprocessing parameters

Allele Frequency Matrix Preprocessing parameters control how the AFM generated after raw reads pre-processing
is filtered and how cell genotypes are assigned.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--filtering` | MT-SNVs filtering method [`null`, `MiTo`, `MQuad`] | `MiTo` |
| `--filter_dbs` | MT-SNVs database filtering | `true` |
| `--spatial_metrics` | Spatial metrics calculation | `false` |
| `--filter_moran` | Spatially segregated MT-SNVs filtering | `true` |
| `--min_cov` | Min MT-SNV site coverage | `5` |
| `--min_var_quality` | Min (average) ALT allele base-calling quality | `30` |
| `--min_frac_negative` | Min fraction of negative cells | `0.2` |
| `--min_n_positive` | Min n of positive cells | `5` |
| `--af_confident_detection` | AF threshold for "confident" detection | `0.02` |
| `--min_n_confidently_detected` | Min n of confidently detected cells | `2` |
| `--min_mean_AD_in_positives` | Min mean AD in +cells | `1.25` |
| `--min_mean_DP_in_positives` | Min mean DP in +cells | `25` |
| `--t_prob` | Probability threshold for MiTo genotyping | `0.7` |
| `--min_AD` | Minimum allelic depth | `2` |
| `--min_cell_prevalence` | Min MT-SNV prevalence for MiTo genotyping | `0.05` |
| `--t_vanilla` | AF threshold for vanilla genotyping | `0` |
| `--bin_method` | Genotyping method [`MiTo`, `vanilla`] | `MiTo` |
| `--min_n_var` | Minimum number of variants per cell | `1` |

`--min_cov`, `--min_var_quality`, `--min_frac_negative`, `--min_n_positive`, `--af_confident_detection`, `--min_n_confidently_detected`, `--min_mean_AD_in_positives`, `--min_mean_DP_in_positives` are active if `--filtering` = `MiTo`. Otherwise all AFM characters are retained (`--filtering` = `null`) or MT-SNVs filtering is performed with the MQuad method (`--filtering` = `MQuad`). See [MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics](https://doi.org/10.1101/2025.06.17.660165) Supplementary Informations for details on MT-SNVs filering and genotyping.

### Phylogeny Reconstruction parameters

Phylogeny Reconstruction parameters control phylogeny reconstruction from filtered character matrices.

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--distance_metric` | Distance metric [`weighted_jaccard`, `jaccard`, `correlation`, `cosine`] | `weighted_jaccard` |
| `--tree_algorithm` | Algorithm [`cassiopeia`, `iqtree`, `mpboot`] | `cassiopeia` |
| `--cassiopeia_solver` | Solver [`UPMGA`, `NJ`, `spectral`, `shared_muts`, `greedy`, `max_cut`] | `UPMGA` |
| `--n_boot_replicates` | Bootstrap replicates | `100` |
| `--boot_strategy` | Bootstrap strategy [`feature_resampling`, `jacknife`] | `feature_resampling` |
| `--frac_char_resampling` | % resampled characters | `0.8` |
| `--support_method` | Support method [`tbe`, `fbp`] | `tbe` |
| `--annotate_tree` | Tree annotation method | `MiTo` |
| `--max_fraction_unassigned` | Tree annotation method | `0.1` |

See [MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics](https://doi.org/10.1101/2025.06.17.660165) Supplementary Informations for details.

## Configuration
Any custom configuration can be passed with the `-c <user.config>` syntax. See the [`config/user.config`](config/user.config) file.

## Examples
Exmples 

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
