# Full parameters

Defaults are set in `nextflow.config` and typed in `nextflow_schema.json`. Give parameters on the command line (`--name value`) or in a JSON file (`-params-file params.json`).

## Input and output

| Parameter | Default | Description |
|---|---|---|
| `--input` | | Samplesheet CSV (`sample,fastq_dir,bam,cell_barcodes`), one library per run. Required. See [Quickstart](quickstart.md). |
| `--assay` | `MAESTER` | `MAESTER`, `mtscATAC` or `ReDeeM`. |
| `--output_folder` | | Folder where the tables are written. Use one folder per library. Required. |

## Reference

| Parameter | Default | Description |
|---|---|---|
| `--masked_ref` | | NUMT-masked genome FASTA, uncompressed. Required. See [Making a masked reference](masked_reference.md). |
| `--string_MT` | `chrM` | Name of the MT contig in the FASTA (and in the BAM). |

## Single-cell tables (all assays)

| Parameter | Default | Description |
|---|---|---|
| `--barcode_tag` | `CB` | BAM tag holding the cell barcode of the MT BAM. |
| `--batch_size` | `auto` | Cells per task. An integer, or `auto`: about 75 tasks per library, of 100 to 1000 cells each. It only sets how the library is divided in tasks; the memory of a task is that of a single cell. |

## MAESTER resources

| Parameter | Default | Description |
|---|---|---|
| `--build_STAR_index` | `true` | Build the STAR index from `--masked_ref` and `--gtf`. |
| `--prebuilt_STAR_index` | | Path to a prebuilt STAR index, used with `--build_STAR_index false`. See [Building a STAR index](star_index.md). |
| `--gtf` | GENCODE v32 | Annotation GTF (path or URL). |
| `--tenx_whitelist` | 10x 3M-february-2018 | 10x barcode whitelist (path or URL). |

## MAESTER (maegatk)

| Parameter | Default | Description |
|---|---|---|
| `--maegatk_umi_tag` | `UB` | UMI tag of the MT BAM. |
| `--maegatk_NHmax` | `2` | Maximum NH (alignments per read). |
| `--maegatk_NMmax` | `15` | Maximum NM (edit distance). |
| `--maegatk_min_reads` | `3` | Reads per UMI family for the fgbio consensus. |
| `--maegatk_base_qual` | `0` | Minimum base quality. |
| `--maegatk_alignment_quality` | `0` | Minimum mapping quality. |
| `--maegatk_proper_paired` | `false` | Keep only properly paired reads. |

## ATAC alignment (mtscATAC, ReDeeM)

| Parameter | Default | Description |
|---|---|---|
| `--atac_min_mapq` | `30` | MAPQ of the MT reads kept after bwa-mem2. |
| `--atac_barcode_rule` | `auto` | Where the 16-nt barcode sits in the barcode read: `auto` (sampled against the cell list), `fwd:<offset>` or `rc:<offset>` (reverse complement). |
| `--barcode_space` | `ATAC` | Space of the barcodes in `cell_barcodes`: `ATAC` (as in the reads) or `RNA` (Cell Ranger ARC, for Multiome/ReDeeM). |
| `--multiome_whitelist_rna` | RedeemV, commit `5c34940` | 10x Multiome RNA whitelist, used with `--barcode_space RNA`. |
| `--multiome_whitelist_atac` | RedeemV, commit `5c34940` | 10x Multiome ATAC whitelist, used with `--barcode_space RNA`. |

## mtscATAC (mgatk)

mgatk runs in `bcall` mode, on one BAM per cell.

| Parameter | Default | Description |
|---|---|---|
| `--mgatk_umi_tag` | empty | UMI tag for duplicate removal. Empty: duplicates are found by coordinates only. |
| `--mgatk_keep_duplicates` | `false` | Keep duplicates. |
| `--mgatk_handle_overlap` | `true` | Handle overlapping mates. |
| `--mgatk_proper_paired` | `true` | Keep only properly paired reads. |
| `--mgatk_base_qual` | `0` | Minimum base quality. |
| `--mgatk_alignment_quality` | `30` | Minimum alignment quality. |
| `--mgatk_NHmax` | `1` | Maximum NH. |
| `--mgatk_NMmax` | `4` | Maximum NM. |

## ReDeeM (RedeemV)

| Parameter | Default | Description |
|---|---|---|
| `--redeem_base_quality` | `30` | Minimum base quality for the molecule consensus. |
| `--redeem_strand_bias_correction` | `storey` | `storey` (Storey q-values, as in RedeemV) or `none` (raw binomial p-values). |
| `--redeem_strand_bias_p` | `0.01` | Threshold on the p-value (`none`) or on the q-value (`storey`). |
| `--redeem_strand_bias_log2fc` | `1` | A biased variant is removed if \|log2(plus/minus)\| is above this value. |

## Generic

| Parameter | Default | Description |
|---|---|---|
| `--help` | `false` | Show the help and exit. |
| `--version` | `false` | Show the version and exit. |
