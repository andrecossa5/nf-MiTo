// bench lineage inference

// Include here
nextflow.enable.dsl = 2
include { MITO_BENCH } from "./modules/mito.nf"
include { CCLONE } from "./modules/cclone.nf"
include { VIREO } from "./modules/vireo.nf"
include { LEIDEN } from "./modules/leiden.nf"


//


// benchmark subworkflow
workflow benchmark {
    
    take:
        ch_jobs 

    main:
    
        MITO_BENCH(ch_jobs)   
        VIREO(ch_jobs)   
        LEIDEN(ch_jobs)   
        CCLONE(ch_jobs)   
        ch_output = MITO_BENCH.out.results.combine(VIREO.out.results, by:[0,1])
            .combine(VIREO.out.results, by:[0,1])
            .combine(LEIDEN.out.results, by:[0,1])
            .combine(CCLONE.out.results, by:[0,1])

    emit:

        pickles = ch_output
        
} 
