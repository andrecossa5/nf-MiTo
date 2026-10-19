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

    OPTIONS:
        --input                  Samplesheet CSV: sample,fastq_dir,bam,cell_barcodes
                                 (exactly one of fastq_dir / bam per row)
        --assay                  MAESTER | mtscATAC | ReDeeM (default: ${params.assay})
        --output_folder          Output directory (required)
        --string_MT              MT contig name (default: ${params.string_MT})
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
