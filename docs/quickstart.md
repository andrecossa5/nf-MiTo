# Quickstart

A run processes **one library** (one row in the samplesheet) and writes the tables straight into `--output_folder`, so give each library its own output folder.

## 1. The samplesheet

`--input` is a CSV with four columns:

```text
sample,fastq_dir,bam,cell_barcodes
donorX,/data/donorX_fastq,,/data/donorX_cells.txt
```

`sample`
: Sample name.

`fastq_dir`
: Folder with the FASTQ files of the library. Give **either** this or `bam`.
  `mtscATAC` and `ReDeeM` expect the cellranger-atac names: `*_R1_*` and `*_R3_*` (genomic reads) and `*_R2_*` (16-nt barcode read). Lanes are concatenated in order.
  `MAESTER` expects `*R1*.fastq.gz` and `*R2*.fastq.gz`.

`bam`
: Coordinate-sorted, indexed BAM with the MT reads and the cell barcode in the `--barcode_tag` tag (default `CB`; MAESTER also needs `UB`). For `ReDeeM` it must come from an alignment to a NUMT-masked genome.

`cell_barcodes`
: Required. A text file with one barcode per line, in the form found in the BAM tag. With FASTQ input, barcodes are matched against the barcode read (up to one mismatch, forward or reverse complement, with or without the `-1` suffix) and reads of other barcodes are dropped. With `--barcode_space RNA` the list is in the RNA space (Cell Ranger ARC) and is converted to ATAC barcodes with the 10x Multiome whitelist pair.

## 2. The parameters

Put the parameters in a JSON file (all of them are described in [Full parameters](parameters.md)):

```json
{
  "assay": "mtscATAC",
  "input": "samplesheet.csv",
  "output_folder": "results/donorX",
  "string_MT": "chrM",
  "masked_ref": "masked.fa"
}
```

## 3. Run

```bash
nextflow run andrecossa5/nf-MiTo -r v2.0.0 -profile singularity -params-file params.json
```

`-profile` is `docker`, `singularity` or `conda`. On a cluster add `-c user.config`, see [Installation](installation.md#cluster-settings). `nextflow run andrecossa5/nf-MiTo -r v2.0.0 --help` lists the parameters.

## Example 1: mtscATAC from a BAM

Download the test BAM bundled with mgatk (37 cells, 50,000 MT reads, 2.7 MB), which also writes the samplesheet:

```bash
bash params/download_mgatk_test_data.sh mgatk_test_data
```

```json
{
  "assay": "mtscATAC",
  "input": "mgatk_test_data/samplesheet.csv",
  "output_folder": "results/cra_v2",
  "string_MT": "chrM",
  "barcode_tag": "CB",
  "masked_ref": "masked.fa"
}
```

With a BAM as input the reads are not aligned: the BAM is split in one pass into batches of cells (`--batch_size`), mgatk (`bcall` mode) counts the cells of each batch one at a time, and the batches are collapsed into the tables.

## Example 2: mtscATAC from FASTQ

```text
sample,fastq_dir,bam,cell_barcodes
GSM9096510,/data/GSM9096510,,/data/GSM9096510_cell_barcodes.txt
```

where `/data/GSM9096510` holds `*_R1_001.fastq.gz`, `*_R2_001.fastq.gz` and `*_R3_001.fastq.gz`. The pipeline concatenates the lanes, reads the cell barcode from R2 (`--atac_barcode_rule auto` picks the offset and orientation by sampling the reads against your cell list), aligns the reads with bwa-mem2 to the masked genome, keeps the MT pairs (`--atac_min_mapq`) and then runs mgatk as in Example 1.

## Example 3: ReDeeM from FASTQ

Download the example library of RedeemV (250,000 read triplets, 10x Multiome ATAC enriched for chrM, and the Cell Ranger ARC cell list), which also writes the samplesheet:

```bash
bash params/download_redeemv_test_data.sh redeemv_test_data
```

The cell list is in the RNA barcode space (Cell Ranger ARC), so we set `--barcode_space RNA`:

```json
{
  "assay": "ReDeeM",
  "input": "redeemv_test_data/samplesheet.csv",
  "output_folder": "results/quickstart",
  "string_MT": "chrM",
  "barcode_space": "RNA",
  "masked_ref": "masked.fa"
}
```

The alignment is the same as for mtscATAC. Then `bin/redeem_consensus.py` builds the molecule consensus genotypes cell by cell (`--redeem_base_quality`), at four stringency levels (`Total`, `VerySensitive`, `Sensitive`, `Specific`), with the rules of RedeemV's `mitoConsensus.py` (commit `5c34940`; output checked to be identical). A **strand-bias filter** then runs per variant, pooled over all cells of the library: a two-sided binomial test of the plus/minus molecule counts against 0.5, converted to Storey q-values (`--redeem_strand_bias_correction storey`, as in RedeemV) or used as raw p-values (`none`). A variant is blacklisted if the p/q-value is below `--redeem_strand_bias_p` and the strands differ by more than `--redeem_strand_bias_log2fc`; its rows are removed unless they have duplex support. Every tested variant is listed in `StrandBiasTable.tsv.gz`.

## Example 4: MAESTER

MAESTER reads are aligned with STARsolo (barcode 1-16, UMI 17-28, `1MM_multi_Nbase_pseudocounts`, `1MM_CR`, `EmptyDrops_CR`: the layout of the MAESTER protocol, fixed), so a STAR index is needed. It is built once from `--masked_ref` and the GENCODE annotation (`--gtf`), or taken from `--prebuilt_STAR_index` with `--build_STAR_index false` (see [Building a STAR index](star_index.md)). `cell_barcodes` lists the barcodes of the cells to analyse, e.g. the cells called on the GEX data.

```json
{
  "assay": "MAESTER",
  "input": "samplesheet.csv",
  "output_folder": "results/donorX",
  "masked_ref": "masked.fa",
  "build_STAR_index": true
}
```

Per cell, maegatk filters the reads (`--maegatk_NHmax`, `--maegatk_NMmax`), builds UMI consensus reads with fgbio (`--maegatk_min_reads`), re-maps them and counts the alleles (`--maegatk_base_qual`, `--maegatk_alignment_quality`).

## Output

```
<output_folder>/
├── A.txt.gz  C.txt.gz  G.txt.gz  T.txt.gz  coverage.txt.gz  depth.txt.gz   # MAESTER, mtscATAC
├── refAllele.txt                                                           # mtscATAC only
├── RawGenotypes.<level>.StrandBalance.gz  QualifiedTotalCts.gz  StrandBias*   # ReDeeM
└── resources/        # reference index, whitelists: reusable by other runs
```

## Load the tables in MiTo

```python
import mito as mt
afm = mt.io.make_afm('results/donorX', 'mtscATAC', 'mgatk')
# ('MAESTER', 'maegatk'), ('ReDeeM', 'redeem-v')
```

## Troubleshooting

Memory errors
: Resources are set per process in `config/base.config`; override them in your `-c` file. Per-cell memory is small and independent of `--batch_size`, which only trades the number of tasks (scheduler load) against the work lost when a task is retried.

Samplesheet errors
: The CSV is checked at launch and the message names the row and column.

Barcodes not found
: `cell_barcodes` must match the BAM tag (or the reads, for FASTQ input). Use `--barcode_space RNA` for Cell Ranger ARC lists. If barcode detection from the reads fails, force the window with `--atac_barcode_rule rc:<offset>` or `fwd:<offset>`.
