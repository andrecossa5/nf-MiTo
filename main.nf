// nf-MiTo
nextflow.enable.dsl = 2

// Preprocess MT
include { preprocess } from "./subworkflows/preprocess/main"
include { createPreprocessingChannel } from "./subworkflows/preprocess/main"

// Preprocess GBC
include { bulk_gbc } from "./subworkflows/bulk_gbc/main"
include { sc_gbc } from "./subworkflows/sc_gbc/main"  

// Phylo inference: tune, explore, phylo
include { afm_preprocess } from "./subworkflows/afm_preprocessing/main"
include { createAFMChannel } from "./subworkflows/afm_preprocessing/main"
include { build_tree } from "./subworkflows/tree_building/main"
include { process_tree } from "./subworkflows/tree_processing/main"
include { hyper_tuning } from "./subworkflows/hyper_tuning/main"
include { explore_spaces } from "./subworkflows/exploring/main"

//

//----------------------------------------------------------------------------//
// nf-MiTo complete entry-points
//----------------------------------------------------------------------------//

// Raw data pre-processing
workflow preprocess_raw {

    ch_preprocessing = createPreprocessingChannel()
    preprocess(ch_preprocessing)

}

//----------------------------------------------------------------------------//

// Lineage inference
workflow tune {

    ch_jobs = createAFMChannel()
    hyper_tuning(ch_jobs)

}

//

workflow explore {

    ch_jobs = createAFMChannel()
    explore_spaces(ch_jobs)

}

//

workflow infer {

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

    ch_preprocessing = createPreprocessingChannel()
    preprocess(ch_preprocessing)
    afm_preprocess(preprocess.out.afm)
    build_tree(afm_preprocess.out.input)
    process_tree(afm_preprocess.out.input, build_tree.out.final_tree)

}


//
