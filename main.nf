// nf-MiTo

nextflow.enable.dsl = 2
include { preprocess } from "./subworkflows/preprocessing/main"
include { build_tree } from "./subworkflows/tree_building/main"
include { process_tree } from "./subworkflows/tree_processing/main"
include { hyper_tuning } from "./subworkflows/hyper_tuning/main"
include { explore_spaces } from "./subworkflows/exploring/main"

//

//----------------------------------------------------------------------------//
// phylo_inference workflow main entry point
//----------------------------------------------------------------------------//

workflow tune {

    ch_jobs = Channel.fromPath(params.path_input)
        .splitCsv(header: true)
        .map { row -> [ row.job_id, row.sample, row.ch_matrix ]}
    hyper_tuning(ch_jobs)
    hyper_tuning.out.summary.view()

}

//

workflow explore {

    ch_jobs = Channel.fromPath(params.path_input)
        .splitCsv(header: true)
        .map { row -> [ row.job_id, row.sample, row.ch_matrix ]}
    explore_spaces(ch_jobs)
    explore_spaces.out.plots.view()

}

//

workflow phylo {

    ch_jobs = Channel.fromPath(params.path_input)
        .splitCsv(header: true)
        .map { row -> [ row.job_id, row.sample, row.ch_matrix ]}
    preprocess(ch_jobs)
    build_tree(preprocess.out.input)
    process_tree(preprocess.out.input, build_tree.out.final_tree)
    process_tree.out.metrics.view()

}

//

// Mock

// Default message
workflow {
    
    println "\n"
    println "Hi there! This is the new version of the phylo_inference pipeline from the MiTo toolkit."
    println "Usage: nextflow run main.nf -c <config> -params-file <params> -profile <profile> -entry <entry>"
    println "See https://github.com/.../main.nf ./config and ./params for configurations and options available."
    println "N.B. This is a BETA version under active development."
    println "\n"

}
