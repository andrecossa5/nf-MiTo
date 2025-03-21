// nf-MiTo
nextflow.enable.dsl = 2

// Preprocess 10x and MT
include { preprocess } from "./subworkflows/preprocess/main"
include { createPreprocessingChannel } from "./subworkflows/preprocess/main"

// Phylo inference: tune, explore, phylo
include { afm_preprocess } from "./subworkflows/afm_preprocess/main"
include { createAFMChannel } from "./subworkflows/afm_preprocess/main"
include { build_tree } from "./subworkflows/tree_build/main"
include { process_tree } from "./subworkflows/tree_process/main"
include { tune } from "./subworkflows/tune/main"
include { explore } from "./subworkflows/explore/main"

//

//----------------------------------------------------------------------------//
// nf-MiTo complete entry-points
//----------------------------------------------------------------------------//

// Raw data pre-processing
workflow PREPROCESS {

    ch_preprocessing = createPreprocessingChannel()
    preprocess(ch_preprocessing)

}

//----------------------------------------------------------------------------//

// Lineage inference
workflow TUNE {

    ch_jobs = createAFMChannel()
    tune(ch_jobs)

}

//

workflow EXPLORE {

    ch_jobs = createAFMChannel()
    explore(ch_jobs)

}

//

workflow INFER {

    ch_jobs = createAFMChannel()
    afm_preprocess(ch_jobs)
    build_tree(afm_preprocess.out.input)
    process_tree(afm_preprocess.out.input, build_tree.out.final_tree)

}


//


//----------------------------------------------------------------------------//
// nf-MiTo main workflow
//----------------------------------------------------------------------------//

workflow {
    
    println "\n"
    println "This is nf-MiTo, the Nextflow pipeline for MT-SNVs-based scLT."
    println "Usage: nextflow run main.nf -c <config> -params-file <params> -profile <profile> -entry <entry>"
    println "See https://github.com/andrecossa5/nf-MiTo for all configurations and options available."
    println "\n"

    ch_pp = createPreprocessingChannel()
    preprocess(ch_pp)
    afm_preprocess(preprocess.out.afm)
    build_tree(afm_preprocess.out.input)
    process_tree(afm_preprocess.out.input, build_tree.out.final_tree)

}


//
