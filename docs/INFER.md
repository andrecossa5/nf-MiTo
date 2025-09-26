# INFER tutorial

nf-MiTo implements an automated workflow (INFER entrypoint, `-entry INFER` option) for lineage inference (cell phylogeny and clones) from pre-processed Allele Frequency Matrices (AFMs). Here, we will demonstrate how to perform such inferences using data from the `MDA_clones` sample, a mixture of lentivirally-barcoded single-cell derived colonies, profiled with a modified version of the MAESTER protocol to include these additional genetic barcodes as ground truth lineage labels (see our [pre-print](https://doi.org/10.1101/2025.06.17.660165) for more details on the experimental design and sequencing strategy).

The whole workflow should take ~... minutes.

# Prerequisites

The clonal inference workflow is not particularly memory intensive (i.e., a modern laptop with 16GB RAM and 8 cpus can handle all processes). However, being on a HPC computing cluster can significantly speed up all parallel operations (e.g., character matrices bootstrapping, distance calculations and tree inference), with multiple processes executed simoultaneously. 
We would also need one between one between docker/apptainer/singularity in our $PATH, together with Nextflow.
We can check this with:

```bash
apptainer 
nextflow
```

**N.B.**: In this tutorial, we will align reads to a genome index built from scratch. nf-MiTo will fetch a genome reference and its annotations from Gencode, and build a STAR index following best practices for single-cell and MT-reads pre-processing (see in [10x documentation](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-3p-references) and in [Lareau et al.](10.1038/s41596-022-00795-3)). This genome index will be available for further runs (through the <span style="white-space: nowrap;">`--prebuilt_STAR_index`</span> parameter). Alternatively, one can:

1. Specify for a different genome and matched annotations (i.e., <span style="white-space: nowrap;">`--reference_genome`</span>, <span style="white-space: nowrap;">`--gtf`</span> parameters)
2. Pre-build a custom STAR index (see [`create_STAR_index.md`](create_STAR_index.md) tutorial).

