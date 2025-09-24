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

Moreover, nf-MiTo supports lineage inference (step 2-5) from pre-processed character matrices of other lineage tracing systems, representing a novel unifying framework for scLT data analysis.

See our recent pre-print [MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics](https://doi.org/10.1101/2025.06.17.660165) for detailed explanations and benchmarks.

### Key Features

- 📁 **Multiple raw data input formats**: raw FASTQ files mitochondrial BAM files
- 🧬 **Multiple lineage tracing systems**: MAESTER, RedeeM, Cas9-based, and scWGS systems
- ⚡  **Highly optimized methods for MT-scLT data analysis**: [MiTo](https://github.com/andrecossa5/MiTo) python package functionalities, at scale
- 🌳 **Flexible phylogeny reconstruction**: different phylogeny reconstruction and boostrap algorithms
- 🎯 **Estensive parameter optimization**: built-in parameter tuning workflow
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

- Run configuration (see ...): `-c <user.config>`, 
- Parameters (see ...): `-params-file <params.json>`, 
- Profiles (see ...): `-profile <chosen_profiles>`

and opt for the main pipeline workflow (end-to-end), or one of the 4 alternative entrypoints (`PREPROCESS`,
`TUNE`, `EXPLORE`, `INFER`, `-entry <chosen_entrypoint>` option).

## Parameters

nf-MiTo implements n=56 parameters controlling one or more of the available entrypoints.
These parameters are grouped in 8 distinct groups. See also `nextflow.config` and `nextflow_schema.json` 
for additional type information. 

### Input/Output parameters

Input/Output Options control nf-MiTo pipeline I/O operations.

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--raw_data_input_type` | Input data type [`fastq`, `fastq, MAESTER`, `mitobam`] | `mitobam` | |
| `--raw_data_input` | Raw data CSV file, preprocessing workflows | `–` | ✅ |
| `--afm_input` | AFM data CSV file, inference workflows | `–` | ✅ |
| `--output_folder` | Output directory path | `–` | ✅ |
| `--path_meta` | Cell metadata file path | `null` | |
| `--path_tuning` | Tuning results file from TUNE workflow | `null` | |

For the default end-to-end workflow or the PREPROCESS entypoint (i.e., `-entry PREPROCESS`), the `--raw_data_input` parameter is required. This parameter points to a CSV file sheet storing IDs and paths
to raw sequencing data. 3 alternative inputs are supported here: 

1. `--raw_data_input_type` = `fastq`. The user need to pre-process both GEX and MT sequening data. MT and GEX FASTQs are passed (`--raw_data_input` parameter) with the following a sample sheet (CSV file):

| sample | fastq_folder | library |
|-----------|-------------|---------|
| `sample_name` | `FASTQs folder path` | `MT` |
| `sample_name` | `FASTQs folder path` | `TENX` |

Each sample is linked to a `fastq_folder` for its `MT` and `TENX` (i.e., GEX) library, respectively. 

2. `--raw_data_input_type` = `fastq, MAESTER`. The user is only interested into MT data pre-processing (i.e., GEX data has been analyzed independently). MT FASTQs are passed (`--raw_data_input` parameter) with the following sample sheet (CSV file):

| sample | fastq_folder | cell_barcodes |
|-----------|-------------|---------|
| `sample_name` | `MT FASTQs folder path` | `cell_barcodes path` |

Each sample is linked to a `fastq_folder` for its `MT` library, and a `cell_barcodes.txt` file with cell barcodes of interests (e.g., qualified cell barcodes from independent GEX data analysis). Only reads associated to these cell barcodes will be filtered during pre-processing of `MT` library.

3. `--raw_data_input_type` = `mitobam`. The user already has aligned MT libraries to the genome reference. These aligned reads (i.e., .bam file) are passed (`--raw_data_input` parameter) with the following a sample sheet (CSV file):

| sample | bam | cell_barcodes |
|-----------|-------------|---------|
| `sample_name` | `MT bam path` | `cell_barcodes path` |

Each sample is linked to a `MT bam file` for its `MT` library, and a `cell_barcodes.txt` file with cell barcodes of interests.

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

### Reference Genome parameters

Input/Output Options control nf-MiTo pipeline I/O operations.

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--ref` | Reference genome directory | - | ✅ |
| `--string_MT` | Mitochondrial chromosome identifier | `chrM` | |
| `--whitelist` | 10x v3 whitelist file | - | ✅ |

### scLT System parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--scLT_system` | System type [MAESTER, RedeeM, Cas9, scWGS] | `MAESTER` | |
| `--pp_method` | Preprocessing method | `maegatk` | |

### Sequencing Data Preprocessing parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--CBs_chunk_size` | Cell barcode chunk size for processing | `3000` | |
| `--fgbio_min_reads_mito` | Minimum reads required for mitochondrial consensus | `3` | |
| `--fgbio_base_quality` | Minimum base quality for consensus calling | `30` | |

### Cell Filtering parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--min_nUMIs` | Minimum UMIs per cell | `500` | |
| `--min_n_genes` | Minimum genes per cell | `250` | |
| `--max_perc_mt` | Maximum mitochondrial read percentage | `0.15` | |
| `--min_cell_number` | Minimum cells with variant | `5` | |
| `--min_cov` | Minimum coverage | `5` | |
| `--min_var_quality` | Minimum variant quality | `30` | |

### Allele Frequency Matrix Preprocessing parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--min_n_positive` | Minimum positive cells | `5` | |
| `--af_confident_detection` | AF threshold for confident detection | `0.02` | |
| `--min_n_confidently_detected` | Minimum confidently detected cells | `2` | |
| `--min_mean_AD_in_positives` | Minimum mean allelic depth | `1.25` | |
| `--t_prob` | Probability threshold | `0.7` | |
| `--min_AD` | Minimum allelic depth | `2` | |
| `--min_cell_prevalence` | Minimum cell prevalence | `0.05` | |
| `--bin_method` | Binarization method | `MiTo` | |

### Phylogeny Reconstruction parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--distance_metric` | Distance metric [weighted_jaccard, jaccard, hamming, cosine] | `weighted_jaccard` | |
| `--tree_algorithm` | Algorithm [cassiopeia, iqtree, mpboot] | `cassiopeia` | |
| `--cassiopeia_solver` | Solver [UPMGA, NJ, UPGMA] | `UPMGA` | |
| `--n_boot_replicates` | Bootstrap replicates | `100` | |
| `--boot_strategy` | Bootstrap strategy [feature_resampling, cell_resampling] | `feature_resampling` | |
| `--lineage_column` | Metadata column for lineage annotation | `null` | |
| `--annotate_tree` | Tree annotation method | `MiTo` | |

### Benchmarking parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--maxK` | Maximum number of clusters for vireo | `15` | |


## Configuration





























## Entry points, and associated outputs

Considering 
nf-MiTo provides several entry points covering different analysis scenarios:


### 1 **Default** - end-to-end data workflow
```bash
nextflow run main.nf  \

```

### 2. **PREPROCESS** - raw sequencing data preprocessing only
Processes raw sequencing data through quality control, alignment, and variant calling to generate allele frequency matrices.

```bash
nextflow run main.nf -entry PREPROCESS \

```

### 3. **TUNE** - critical parameter optimization
Systematically tests different parameter combinations to optimize variant detection for your specific dataset.

```bash
nextflow run main.nf -entry TUNE \

```

### 4. **EXPLORE** - alternative MT-SNVs spaces visualization
Generates comprehensive visualizations and quality metrics for mitochondrial variant space exploration.

```bash
nextflow run main.nf -entry EXPLORE \

```

### 5. **INFER** - lineage inference
Performs phylogenetic reconstruction and tree annotation using pre-computed or optimized parameters.

```bash
nextflow run main.nf -entry INFER \

```

## Input Data Formats

### Raw Data Input (`--raw_data_input`)

For preprocessing workflows, provide a CSV file with the following structure:

**For `--raw_data_input_type fastq`:**
```csv
sample,r1,r2
sample1,/path/to/sample1_R1.fastq.gz,/path/to/sample1_R2.fastq.gz
sample2,/path/to/sample2_R1.fastq.gz,/path/to/sample2_R2.fastq.gz
```

**For `--raw_data_input_type mitobam`:**
```csv
sample,path
sample1,/path/to/sample1_mito.bam
sample2,/path/to/sample2_mito.bam
```

### AFM Input (`--afm_input`)

For direct analysis workflows (TUNE, INFER, EXPLORE, BENCH), provide a CSV file:

```csv
job_id,sample,afm
job1,sample1,/path/to/sample1_afm.csv
job2,sample2,/path/to/sample2_afm.csv
```









## Parameter Tuning Strategy

The TUNE workflow systematically tests parameter combinations to optimize variant detection:

### Tunable Parameters
When using the TUNE entry point, specify arrays for key parameters:

```bash
# Example parameter file for tuning
{
    "min_n_positive": [3, 5, 7],
    "af_confident_detection": [0.01, 0.02, 0.03],
    "t_prob": [0.6, 0.7, 0.8],
    "bin_method": ["MiTo", "vanilla"]
}
```

### Tuning Workflow
1. **Grid Search**: Tests all parameter combinations
2. **Quality Metrics**: Evaluates variant detection quality
3. **Optimal Selection**: Identifies best parameter set
4. **Results Export**: Provides recommendations for INFER

## Examples

### Example 1: Complete Pipeline from FASTQ

```bash
nextflow run main.nf \
    -profile docker \
    --raw_data_input_type fastq \
    --raw_data_input samples.csv \
    --ref /data/reference/hg38 \
    --whitelist /data/barcodes/3M-february-2018.txt \
    --scLT_system MAESTER \
    --output_folder complete_analysis
```

### Example 2: Parameter Tuning for Optimal Detection

```bash
nextflow run main.nf \
    -entry TUNE \
    -profile docker \
    -params-file tune_params.json \
    --afm_input afm_samples.csv \
    --output_folder parameter_tuning
```

**tune_params.json:**
```json
{
    "afm_input": "afm_samples.csv",
    "output_folder": "parameter_tuning",
    "min_n_positive": [3, 5, 7, 10],
    "af_confident_detection": [0.01, 0.02, 0.03],
    "min_n_confidently_detected": [2, 3, 4],
    "t_prob": [0.6, 0.7, 0.8],
    "bin_method": ["MiTo", "vanilla"]
}
```

### Example 3: Inference with Custom Parameters

```bash
nextflow run main.nf \
    -entry INFER \
    -profile docker \
    --afm_input afm_samples.csv \
    --path_tuning /results/tuning/optimal_params.csv \
    --lineage_column cell_type \
    --distance_metric weighted_jaccard \
    --tree_algorithm cassiopeia \
    --n_boot_replicates 1000 \
    --output_folder lineage_inference
```

### Example 4: Benchmarking Different Methods

```bash
nextflow run main.nf \
    -entry BENCH \
    -profile docker \
    --afm_input afm_samples.csv \
    --maxK 20 \
    --output_folder benchmarking
```

## Output Structure

```
output_folder/
├── preprocessing/          # Raw data processing results
│   ├── qc/                # Quality control metrics
│   ├── alignment/         # BAM files and indices
│   └── variants/          # Variant calling results
├── afm_preprocessing/     # AFM generation and filtering
│   ├── matrices/          # Raw and filtered AFMs
│   ├── metrics/           # Quality metrics
│   └── visualization/     # QC plots
├── tuning/               # Parameter optimization results
│   ├── grid_search/      # All parameter combinations
│   ├── metrics/          # Performance metrics
│   └── optimal/          # Best parameter sets
├── phylogeny/            # Tree reconstruction
│   ├── trees/            # Phylogenetic trees
│   ├── bootstrap/        # Bootstrap support
│   └── annotation/       # Annotated trees
├── benchmarking/         # Method comparison
│   ├── clustering/       # Different clustering results
│   ├── metrics/          # Benchmark metrics
│   └── comparison/       # Comparative analysis
└── reports/              # Summary reports and plots
    ├── html/             # Interactive reports
    ├── figures/          # Publication-ready figures
    └── tables/           # Summary statistics
```

## Best Practices

### 1. **Start with Parameter Tuning**
For new datasets, always run TUNE first to optimize variant detection parameters.

### 2. **Quality Control**
- Monitor cell filtering metrics
- Check variant detection quality
- Validate phylogenetic reconstruction

### 3. **Resource Management**
- Use appropriate computing resources for dataset size
- Consider chunking large datasets
- Monitor memory usage for tree reconstruction

### 4. **Reproducibility**
- Use versioned containers
- Document parameter choices
- Save configuration files

## Execution Profiles

### Docker (Recommended)
```bash
nextflow run main.nf -profile docker
```

### Singularity (HPC environments)
```bash
nextflow run main.nf -profile singularity
```

### Conda (Local development)
```bash
nextflow run main.nf -profile conda
```

### Local (No containers)
```bash
nextflow run main.nf -profile local
```

## Troubleshooting

### Common Issues

1. **Memory Errors**: Increase memory allocation in configuration
2. **Missing Dependencies**: Ensure correct execution profile
3. **Input Format Errors**: Validate CSV file structure
4. **Reference Genome**: Verify genome index completeness

### Support

- **Documentation**: [GitHub Repository](https://github.com/andrecossa5/nf-MiTo)
- **Issues**: Report bugs via GitHub Issues
- **Discussions**: Community support via GitHub Discussions

## Citation

If you use nf-MiTo in your research, please cite:

```
[Citation information to be added]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
