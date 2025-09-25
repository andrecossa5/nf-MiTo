# PREPROCESS tutorial


# Prerequisites

First, we make sure that we are in a conda/mamba environment with docker/apptainer/singularity and Nextflow installed.

```bash
apptainer 
nextflow
```

**N.B.**: This tutorial assumes that we want to use nf-MiTo built-in genome reference, transcripts annotation, and 10x cellular
barcodes whitelist. If this is not be the case (e.g., a specific reference genome is needed), we would need to built our own STAR index (see [`create_STAR_index.md`](create_STAR_index.md) tutorial).

# Download data

We start by downloding test data from zenodo

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

We have downloaded paired-end FASTQ files for the MAESTER (MT) and TENX (GENX) libraries of n=5 test cells. These will be our test data for the PREPROCESS entrypoint.

# Prep input file

We prep our `--raw_data_input` CSV file:

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

Here, we will *not* apply any custom `user.config`. We will just apply the `docker` and `local` profiles (defined in [`config/base.config`](config/base.config)). With this configs, each process will be run by the local executor in its dedicated docker container (pulled locally from [Dockerhub](https://hub.docker.com/)). 

However, we can pass custom `user.configs` with the `-c` option, and make our run compatible with virtually any computing environment. For example, check [`config/user.config`](config/user.config). This custom config file specify for additional process directives profiles to run nf-MiTo on an HPC cluster with singularity and the SLURM job scheduler. For furtherly customized options, see [profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles).

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

We are now ready to launch `nf-MiTo` PREPROCESS. 

```bash
nextflow run main.nf -profile docker,local -params-file user_params.json -entry PREPROCESS
```

# Inspect output

The output folder `<path to output folder>` should contain all the following outputs:

```
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

...


# Alternative --raw_data_input_type options






















