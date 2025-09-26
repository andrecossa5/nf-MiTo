# INFER tutorial

nf-MiTo implements an automated workflow (INFER entrypoint, `-entry INFER` option) for lineage inference (cell phylogeny and clones) from pre-processed Allele Frequency Matrices (AFMs. See the [PREPROCESSING tutorial](./PREPROCESS.md) for raw sequencing data pre-processing). Here, we will demonstrate how to perform such inferences using data from the `MDA_clones` sample, a mixture of lentivirally-barcoded single-cell derived colonies, profiled with a modified version of the MAESTER protocol to include these additional genetic barcodes as ground truth lineage labels (see our [pre-print](https://doi.org/10.1101/2025.06.17.660165) for more details on the experimental design and sequencing strategy).

The whole workflow should take ~15 minutes (on an HPC environment).

# Prerequisites

The clonal inference workflow is not particularly memory intensive (i.e., a modern laptop with 16GB RAM and 8 cpus can handle all processes). However, being on a HPC computing cluster can significantly speed up all parallel operations (e.g., character matrices bootstrapping, distance calculations and tree inference), with multiple processes executed simoultaneously. 
As mentioned in the other tutorials, our only requirements include one between docker/apptainer/singularity in our $PATH, together with Nextflow.
We can check their availability with:

```bash
apptainer 
nextflow
```

# Download data

We start by downloding test data from zenodo. This will include the AFM (and corresponding cell metadata)
output of the PREPROCESS entrypoint (default options):

```bash
wget https://zenodo.org/records/17209529/files/data_test.tar.gz
```

# Prep data

Untar test_data:

```bash
tar -xzf data_test.tar.gz
```

Our `data_test` folder should look like:

```
data_test
├── afm_unfiltered.h5ad
└── cells_meta.csv
```

We have downloaded our test dataset. We are now ready to prep our `afm_input` file.

# Prep input file

We prep our <span style="white-space: nowrap;">`--afm_input`</span> CSV file following this template:

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
<td align="left"><code>test_0.1</code></td>
<td align="left">&lt;MDA_clones&gt;</td>
<td align="left"><code>path to AFM</code></td>
</tr>
</tbody>
</table>

# Set parameters

We prep a `user_params.json` file with our custom parameters:

```json
{  
    "afm_input" : "<path to --afm_input .csv file>",
    "output_folder" : "<path to output folder>",
    "path_meta" : "<path to cells metadata .csv file>"     
}
```

# Set configs

See [PREPROCESS tutorial](PREPROCESS.md) for config settings.


# Launch pipeline

We can now launch `nf-MiTo` INFER:

```bash
nextflow run main.nf -profile docker,local -params-file user_params.json -entry INFER
```

# Inspect output

After a successfull run, the output folder `<path to output folder>` should contains all the following outputs:

```
├── INFER
│   └── MDA_clones
│       └── test_0.1
│           ├── afm_filtered.h5ad
│           ├── annotated_tree.pickle
│           └── tree_metrics.csv
```

Each sample in the run (`MDA_clones`, in this case) has its own folder, with as many subfolders (named after their job_id) as the AFMs associated to that sample in the `--afm_input` file. This design allows one to process multiple unfiltered character matrices per sample (e.g. different pre-processing pipelines, if necessary), tagging these alternative outputs with a `job_id` tag. For each of these lineage inference jobs, 3 outputs are produced:

* `afm_filtered.h5ad`: the AFM, after cell and MT-SNVs filtering, genotyping and distance calculations.
* `annotated_tree.pickle`: the cell phylogenies, with annotated metadata and discrete clonal labels (if `--annotate_tree` = `MiTo`)
* `tree_metrics.csv`: metrics evaluating several properties of the inferred phylogenies (see our [pre-print](https://doi.org/10.1101/2025.06.17.660165) for details).

# Alternative pre-processing workflows

See all the available options to tune lineage inference in the main [README](../README.md).

