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
            --raw_data_input samples.csv --raw_data_input_type <type> \\
            --output_folder results

    OPTIONS:
        --raw_data_input_type    "fastq, MAESTER" | mitobam   (default: ${params.raw_data_input_type})
        --raw_data_input         CSV with the samples
        --output_folder          Output directory (required)
        --scLT_system            ${params.scLT_system}
        --pp_method              ${params.pp_method}
        --string_MT              MT contig name (default: ${params.string_MT})

    Work in progress: mgatk (mtscATAC) and RedeemV (ReDeeM) are being added.
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
