# nf-MiTo
**Integrated Nextflow pipeline for mitochondrial SNV-based single-cell lineage tracing and multi-omics**

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A522.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-enabled-blue.svg)](https://www.docker.com/)

## Overview

nf-MiTo is a flexible Nextflow pipeline for mitochondrial SNV-based single-cell lineage tracing (MT-scLT). 
The pipeline is specifically focused on [MAESTER](10.1038/s41587-022-01210-8) data, where MT-SNVs are profiled from full lenght 10x cDNA. In particular, nf-MiTo implements end-to-end workflows for coupled analysis of Gene Expression (GEX) and MT-SNVs (MT) data raw sequencing data (.fastq or .bam format):

1. Preprocessing (raw sequencing data)
2. Allele Frequency Matrix (AFM) filtering
3. Cell Genotyping
4. Calculation of cell-cell distances
5. Lineage inference

Moreover, nf-MiTo supports lineage inference (step 2-5) from pre-processed character matrices of other lineage tracing systems, representing a novel unifying framework for scLT data analysis.

### Key Features

- **Multiple raw data input formats**: raw FASTQ files mitochondrial BAM files
- **Multiple lineage tracing systems**: MAESTER, RedeeM, Cas9-based, and scWGS systems
- **Highly optimized methods for MT-scLT data analysis**: MiTo python package functionalities, at scale
- **Flexible phylogeny reconstruction**: different phylogeny reconstruction algorithms, Transfer Bootstrap support
- **Estensive parameter optimization**: built-in parameter tuning workflow
- **Scalable execution**: in local, HPC, and cloud environments

## Quick Start

### Prerequisites

- [Nextflow](https://www.nextflow.io/) (≥22.04.0)
- [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) or [Conda](https://conda.io/)

### Basic Usage

```bash
# Run the full pipeline with mitobam input
nextflow run main.nf \
    -profile docker \
    --raw_data_input samples.csv \
    --output_folder results \
    --ref /path/to/reference

# Run parameter tuning
nextflow run main.nf \
    -entry TUNE \
    -profile docker \
    --afm_input afm_jobs.csv \
    --output_folder tune_results

# Run inference only  
nextflow run main.nf \
    -entry INFER \
    -profile docker \
    --afm_input afm_jobs.csv \
    --output_folder infer_results
```

## Workflow Modes

nf-MiTo provides several entry points for different analysis scenarios:

### 1. **PREPROCESS** - Data Preprocessing Only
Processes raw sequencing data through quality control, alignment, and variant calling to generate allele frequency matrices.

```bash
nextflow run main.nf -entry PREPROCESS \
    --raw_data_input samples.csv \
    --output_folder preprocess_results
```

### 2. **TUNE** - Parameter Optimization
Systematically tests different parameter combinations to optimize variant detection for your specific dataset.

```bash
nextflow run main.nf -entry TUNE \
    --afm_input afm_jobs.csv \
    --output_folder tune_results
```

### 3. **INFER** - Lineage Inference (Default)
Performs phylogenetic reconstruction and tree annotation using pre-computed or optimized parameters.

```bash
nextflow run main.nf -entry INFER \
    --afm_input afm_jobs.csv \
    --output_folder infer_results
```

### 4. **EXPLORE** - Exploratory Analysis
Generates comprehensive visualizations and quality metrics for mitochondrial variant space exploration.

```bash
nextflow run main.nf -entry EXPLORE \
    --afm_input afm_jobs.csv \
    --output_folder explore_results
```

### 5. **BENCH** - Benchmarking
Compares different clustering methods and evaluates pipeline performance across parameter sets.

```bash
nextflow run main.nf -entry BENCH \
    --afm_input afm_jobs.csv \
    --output_folder bench_results
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

## Supported Lineage Tracing Systems

### MAESTER (Default)
Mitochondrial RNA-editing-based lineage tracing with UMI-based consensus calling.

### RedeeM  
Mitochondrial base editing system for lineage tracing.

### Cas9
CRISPR/Cas9-based mitochondrial editing for lineage tracking.

### scWGS
Single-cell whole genome sequencing with mitochondrial variant focus.

## Configuration Parameters

### Essential Parameters

| Parameter | Description | Default | Required |
|-----------|-------------|---------|----------|
| `--output_folder` | Output directory path | - | ✅ |
| `--raw_data_input` | Raw data CSV file | - | For PREPROCESS |
| `--afm_input` | AFM data CSV file | - | For INFER/TUNE/EXPLORE/BENCH |
| `--ref` | Reference genome directory | - | For PREPROCESS |

### Input/Output Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--raw_data_input_type` | Input data type [fastq, "fastq, MAESTER", mitobam] | `mitobam` |
| `--path_meta` | Cell metadata file path | `null` |
| `--path_tuning` | Tuning results file path | `null` |

### Reference Genome Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--string_MT` | Mitochondrial chromosome identifier | `chrM` |
| `--whitelist` | Cell barcode whitelist file | - |

### Lineage Tracing System

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--scLT_system` | System type [MAESTER, RedeeM, Cas9, scWGS] | `MAESTER` |
| `--pp_method` | Preprocessing method | `maegatk` |

### Quality Control & Filtering

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--min_nUMIs` | Minimum UMIs per cell | `500` |
| `--min_n_genes` | Minimum genes per cell | `250` |
| `--max_perc_mt` | Maximum mitochondrial read percentage | `0.15` |
| `--min_cell_number` | Minimum cells with variant | `5` |
| `--min_cov` | Minimum coverage | `5` |
| `--min_var_quality` | Minimum variant quality | `30` |

### Variant Detection (Key Parameters for Tuning)

| Parameter | Description | Default | Tunable |
|-----------|-------------|---------|---------|
| `--min_n_positive` | Minimum positive cells | `5` | ✅ |
| `--af_confident_detection` | AF threshold for confident detection | `0.02` | ✅ |
| `--min_n_confidently_detected` | Minimum confidently detected cells | `2` | ✅ |
| `--min_mean_AD_in_positives` | Minimum mean allelic depth | `1.25` | ✅ |
| `--t_prob` | Probability threshold | `0.7` | ✅ |
| `--min_AD` | Minimum allelic depth | `2` | ✅ |
| `--min_cell_prevalence` | Minimum cell prevalence | `0.05` | ✅ |
| `--bin_method` | Binarization method | `MiTo` | ✅ |

### Phylogeny Reconstruction

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--distance_metric` | Distance metric [weighted_jaccard, jaccard, hamming, cosine] | `weighted_jaccard` |
| `--tree_algorithm` | Algorithm [cassiopeia, iqtree, mpboot] | `cassiopeia` |
| `--cassiopeia_solver` | Solver [UPMGA, NJ, UPGMA] | `UPMGA` |
| `--n_boot_replicates` | Bootstrap replicates | `100` |
| `--boot_strategy` | Bootstrap strategy [feature_resampling, cell_resampling] | `feature_resampling` |
| `--lineage_column` | Metadata column for lineage annotation | `null` |
| `--annotate_tree` | Tree annotation method | `MiTo` |

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
