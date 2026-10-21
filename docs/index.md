# nf-MiTo

nf-MiTo is a [Nextflow](https://www.nextflow.io/) pipeline for the **preprocessing** of mitochondrial single-cell lineage tracing (MT-scLT) data. From the sequencing reads (FASTQ) or an aligned BAM of a library, it produces the single-cell allele tables read by [MiTo](https://github.com/andrecossa5/MiTo) (`mt.io.make_afm`). Filtering, genotyping, trees and clones are done in MiTo, not here.

![nf-MiTo workflow](../image/nf_mito.svg)

Three assays are supported, one tool each:

| `--assay` | Library | Tool (`pp_method` in MiTo) | Alignment |
|---|---|---|---|
| `MAESTER` | 10x GEX + MAESTER MT amplicons | [maegatk](https://github.com/petervangalen/MAESTER-2021) (`maegatk`) | STARsolo |
| `mtscATAC` | 10x scATAC | [mgatk](https://github.com/caleblareau/mgatk) (`mgatk`) | bwa-mem2 |
| `ReDeeM` | 10x Multiome ATAC, ReDeeM protocol | [RedeemV](https://github.com/chenweng1991/redeemV) (`redeem-v`) | bwa-mem2 |

Every assay follows the same two stages:

1. **Reads to MT BAM** (skipped if you give a BAM): one library-level, coordinate-sorted BAM of the chrM reads, with the cell barcode in a tag.
2. **MT BAM to single-cell tables**: the BAM is split in one pass into batches of cells, the assay's tool processes the cells of each batch one at a time (memory is that of one cell), and the batches are collapsed into the tables.

```{toctree}
:maxdepth: 2
:caption: Contents

installation
quickstart
masked_reference
parameters
star_index
```

## Citation

If you use nf-MiTo, please cite *MiTo: tracing the phenotypic evolution of somatic cell lineages via mitochondrial single-cell multi-omics* (Cossa et al., bioRxiv, [doi:10.1101/2025.06.17.660165](https://doi.org/10.1101/2025.06.17.660165)), and the tool you used: maegatk/MAESTER (Miller et al., Nat. Biotechnol. 2022), mgatk (Lareau et al., Nat. Biotechnol. 2021) or RedeemV/ReDeeM (Weng et al., Nature 2024).

nf-MiTo is released under the MIT license.
