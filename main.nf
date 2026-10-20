// nf-MiTo 2.0.0: reads (FASTQ or BAM) -> allelic tables for MiTo
nextflow.enable.dsl = 2

include { preprocess } from "./subworkflows/preprocess/main"
include { createPreprocessingChannel } from "./subworkflows/preprocess/main"


//


if (params.help) {

    println """\
    ==========================================================================
    nf-MiTo 2.0.0: pre-processing of MT-scLT data (reads -> allelic tables)
    ==========================================================================

    USAGE:
        nextflow run main.nf -profile <docker|singularity|conda> \\
            --input samplesheet.csv --assay MAESTER --output_folder results

    INPUT/OUTPUT:
        --input                  Samplesheet CSV: sample,fastq_dir,bam,cell_barcodes
                                 (exactly one of fastq_dir / bam per row, one library per run)
        --assay                  MAESTER | mtscATAC | ReDeeM (default: ${params.assay})
        --output_folder          Output directory (required)

    REFERENCE:
        --masked_ref             NUMT-masked genome FASTA, uncompressed (required)
        --string_MT              MT contig name (default: ${params.string_MT})

    SINGLE-CELL TABLES (all assays):
        --barcode_tag            Cell barcode tag of the MT BAM (default: ${params.barcode_tag})
        --batch_size             Cells per task: integer or auto (default: ${params.batch_size})

    MAESTER:
        --build_STAR_index, --prebuilt_STAR_index, --gtf, --tenx_whitelist   STARsolo resources
        --maegatk_umi_tag, --maegatk_NHmax, --maegatk_NMmax, --maegatk_min_reads,
        --maegatk_base_qual, --maegatk_alignment_quality, --maegatk_proper_paired   maegatk

    mtscATAC / ReDeeM alignment:
        --atac_min_mapq, --atac_barcode_rule, --barcode_space,
        --multiome_whitelist_rna, --multiome_whitelist_atac

    mtscATAC (mgatk):
        --mgatk_umi_tag, --mgatk_keep_duplicates, --mgatk_handle_overlap, --mgatk_proper_paired,
        --mgatk_base_qual, --mgatk_alignment_quality, --mgatk_NHmax, --mgatk_NMmax

    ReDeeM (RedeemV):
        --redeem_base_quality, --redeem_strand_bias_correction, --redeem_strand_bias_p,
        --redeem_strand_bias_log2fc

    Defaults and descriptions: nextflow.config, nextflow_schema.json, README.md
    ==========================================================================
    """.stripIndent()
    exit(0)

}

if (params.version) {
    println "nf-MiTo version: ${workflow.manifest.version}"
    exit(0)
}

if (!params.output_folder) {
    error "Error: --output_folder is required. Use --help for usage information."
}


//


workflow PREPROCESS {

    ch_preprocessing = createPreprocessingChannel()
    preprocess(ch_preprocessing)

}

workflow {

    PREPROCESS()

}
