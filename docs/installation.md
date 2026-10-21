# Installation

## Requirements

- [Nextflow](https://www.nextflow.io/) 24.04.0 or later.
- A container engine (Docker, or Singularity/Apptainer on a cluster) **or** Conda/Mamba.
- An x86_64 machine: bwa-mem2 (`mtscATAC`, `ReDeeM`) and mgatk/maegatk are not built for ARM.
- A NUMT-masked genome FASTA, see [Making a masked reference](masked_reference.md).
- Memory: building the bwa-mem2 index of a human genome needs about 80 GB, once. Per-cell tasks need little.

## Nextflow

```bash
mamba create -n nextflow -c bioconda -c conda-forge nextflow
mamba activate nextflow
nextflow -version
```

## The pipeline

Nextflow pulls the pipeline from GitHub, so there is nothing to install:

```bash
nextflow run andrecossa5/nf-MiTo -r v2.0.0 --help
```

To work on a local copy, clone it and run `nextflow run <path to nf-MiTo>/main.nf`:

```bash
git clone https://github.com/andrecossa5/nf-MiTo.git
```

## Software

nf-MiTo uses one image per tool, built from the recipes in [`envs/`](https://github.com/andrecossa5/nf-MiTo/tree/master/envs). Choose a profile with `-profile`:

| Profile | What it does |
|---|---|
| `docker` | Pulls the images below. |
| `singularity` | Same, converted by Singularity/Apptainer. |
| `conda` | Builds the environments from `envs/*.yml`, or uses pre-installed ones (see below). |

| Label | Image | Used for |
|---|---|---|
| `maegatk` | `acox1/mito-maegatk:2.0.0` | MAESTER, BAM splitting into batches, STARsolo |
| `mgatk` | `acox1/mito-mgatk:2.0.0` | bwa-mem2 alignment, mgatk (v0.9.0) |
| `redeemv` | `acox1/mito-redeemv:2.0.0` | RedeemV consensus and strand-bias filter |

To build an image yourself, from the repository root:

```bash
docker build --platform linux/x86_64 -f docker/Dockerfile_mgatk -t acox1/mito-mgatk:2.0.0 .
```

## Cluster settings

Site settings (executor, queues, scratch, pre-installed conda environments) go in a Nextflow config passed with `-c`. A template is in [`config/user.config`](https://github.com/andrecossa5/nf-MiTo/blob/master/config/user.config):

```bash
nextflow run andrecossa5/nf-MiTo -r v2.0.0 -profile slurm,singularity -c user.config -params-file params.json
```

Per-process resources are set in `config/base.config`; override them in your `-c` file. For example `BWAMEM2_INDEX` requests 80 GB, and `BWAMEM2_MT` 40 GB.

## Loading the tables

The tables are read by [MiTo](https://github.com/andrecossa5/MiTo), which is a separate package (nf-MiTo 2.0.0 does not use it):

```bash
pip install scmito
```
