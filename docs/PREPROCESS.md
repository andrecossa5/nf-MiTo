# PREPROCESS tutorial

nf-MiTo implements an automated workflow (PREPROCESS entrypoint, `-entry PREPROCESS` option) for pre-processing of raw sequencing data from both MAESTER and 10x libraries. Here we will pre-process a minimal test dataset with paired-end FASTQ files for n=5 cells. The whole workflow should take ~45 minutes, most of which (~35) will be needed for building the genome index needed for sequence alignment. 

# Prerequisites

Raw-reads pre-processing requires at least ~50 GB RAM, and at least 0.1-0.5 TB of scratch space depending on the dataset. Note that, regardless of the minimal dataset in this tutorial (<0.5 GB), STAR still needs ~30GB to build a genome index and load it into memory for alignment. 
We would also need one between one between docker/apptainer/singularity in our $PATH, together with Nextflow.
We can check this with:

```bash
apptainer 
nextflow
```

**N.B.**: In this tutorial, we will align reads to a genome index built from scratch. nf-MiTo will fetch a genome reference and its annotations from Gencode, and build a STAR index following best practices for single-cell and MT-reads pre-processing (see in [10x documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references) and in [Lareau et al.](10.1038/s41596-022-00795-3)). This genome index will be available for further runs (through the `--prebuilt_STAR_index` parameter). Alternatively, one can:

1. Specify for a different genome and matched annotations (i.e., `--reference_genome`, `--gtf` parameters)
2. Pre-build a custom STAR index (see [`create_STAR_index.md`](create_STAR_index.md) tutorial).

# Download data

We start by downloding test data from zenodo:

```bash
wget https://zenodo.org/records/17193549/files/raw_data_test.tar.gz
```

# Prep data

Untar test_data:

```bash
tar -xzf raw_data_test.tar.gz
```

Our `test_data` folder should look like:

```
└── test_data
    ├── MAESTER_small
    │   ├── R1_small.fastq.gz
    │   └── R2_small.fastq.gz
    └── TENX_small
        ├── R1_small.fastq.gz
        └── R2_small.fastq.gz
```

We have downloaded our test dataset. We are now ready to prep our `raw_data_input` file.

# Prep input file

We prep our `--raw_data_input` CSV file following this template:

```
| sample | fastq_folder | library |
|-----------|-------------|---------|
| `test` | <path to MAESTER_small> | `MT` |
| `test` | <path to TENX_small> | `TENX` |
```

# Set parameters

We prep a `user_params.json` file with our custom parameters:

```json
{   
    "raw_data_input_type" : "fastq",
    "raw_data_input" : "<path to --raw_data_input .csv file>",
    "output_folder" : "<path to output folder>",
}
```

# Set configs

Here, we will *not* apply any custom `user.config`. We will just apply the `docker` and `local` profiles (defined in [`config/base.config`](config/base.config)). With this configs, each process will be run by our local executor in its dedicated docker container (pulled locally from [Dockerhub](https://hub.docker.com/)). 

However, we can pass custom `user.configs` with the `-c` option, and make our run compatible with virtually any computing environment (check [`config/user.config`](config/user.config) for an example with process directives and profiles designed to run nf-MiTo on an HPC cluster with singularity, and the SLURM job scheduler). For more on run options and configs, see [profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles).

**N:B** to cache Docker/Singularity .img and .sif files across repeated runs, check for your $ENV variables in `~/.bashrc` (or `~/.zshrc`) files:

```
# e.g., Apptainer cache variables 
export NXF_SINGULARITY_CACHEDIR=<apptainer image cache>
export APPTAINER_CACHEDIR=<apptainer image cache>
export APPTAINER_TMPDIR=<local scratch>

# Other variables
# ...
```

# Launch pipeline

We can now launch `nf-MiTo` PREPROCESS:

```bash
nextflow run main.nf -profile docker,local -params-file user_params.json -entry PREPROCESS
```

# Inspect output

After a successfull run, the output folder `<path to output folder>` should contains all the following outputs:

```
└── PREPROCESS
    ├── resources
    │   └── STAR_index
    │       ├── chrLength.txt
    │       ├── chrNameLength.txt
    │       ├── chrName.txt
    │       ├── chrStart.txt
    │       ├── exonGeTrInfo.tab
    │       ├── exonInfo.tab
    │       ├── geneInfo.tab
    │       ├── Genome
    │       ├── genomeParameters.txt
    │       ├── Log.out
    │       ├── SA
    │       ├── SAindex
    │       ├── sjdbInfo.txt
    │       ├── sjdbList.fromGTF.out.tab
    │       ├── sjdbList.out.tab
    │       └── transcriptInfo.tab
    └── test
        ├── adata.h5ad
        ├── afm_unfiltered.h5ad
        ├── cell_barcodes.txt
        ├── MT.preprocess.ouput
        │   └── tables
        │       ├── A.txt.gz
        │       ├── coverage.txt.gz
        │       ├── C.txt.gz
        │       ├── depth.txt.gz
        │       ├── G.txt.gz
        │       └── T.txt.gz
        └── Solo.output
            ├── Aligned.sortedByCoord.out.bam
            ├── Features.stats
            ├── filtered
            │   ├── barcodes.tsv.gz
            │   ├── features.tsv.gz
            │   └── matrix.mtx.gz
            ├── raw
            │   ├── barcodes.tsv.gz
            │   ├── features.tsv.gz
            │   └── matrix.mtx.gz
            └── Summary.csv
```

The `resources` folder store the STAR index we have built (available for further use).

For all the other samples specified in the first columns of the `--raw_data_input_file` (`test` in this case), the `--output_folder` directory contains:

* `adata.h5ad`: AnnData object containing gene expression count matrix and cell metadata from Quality-Controlled (QC) cells. Includes UMI counts per gene, cell filtering statistics, and quality control metrics
* `afm_unfiltered.h5ad`: AnnData object containing the Allele Frequency Matrix (AFM) with mitochondrial variant data from QC cells. Contains base counts (A, C, G, T) at each MT position for each cell, along with coverage and quality metrics
* `cell_barcodes.txt`: Text file listing cell barcodes that passed quality control filters (without sample name suffixes). These barcodes correspond to cells included in both the gene expression and MT datasets
* `MT.preprocess.ouput/`: Directory containing detailed MT preprocessing results:
  * `tables/`: Base-specific count tables (A.txt.gz, C.txt.gz, G.txt.gz, T.txt.gz) with UMI counts for each base at each MT position per cell
  * `coverage.txt.gz`: Per-position coverage depth across all cells
  * `depth.txt.gz`: Total sequencing depth statistics per cell and position
* `Solo.output/`: STAR Solo alignment results for the gene expression (TENX) library:
  * `Aligned.sortedByCoord.out.bam`: Coordinate-sorted BAM file with aligned reads
  * `filtered/`: Cell Ranger-compatible filtered count matrices (barcodes, features, matrix)
  * `raw/`: Unfiltered count matrices including all detected barcodes
  * `Features.stats` and `Summary.csv`: Alignment and quantification statistics


# Alternative pre-processing workflows

Alternative inputs (`--raw_data_input_type` and `--raw_data_input` parameters) and pre-processing tools/pipelines (`--pp_method`) can be specified in `user_params.json`, and passed to the pipeline with the `-params-file` option.
See the [README](../README.md) for details.





















