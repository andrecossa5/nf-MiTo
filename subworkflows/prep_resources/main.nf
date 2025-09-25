// prep_resources workflow

// Include here
nextflow.enable.dsl = 2
include { FETCH_WHITELIST } from "./modules/fetch_whitelist.nf"
include { PREP_GTF } from "./modules/prep_gtf.nf"
include { PREP_GENOME } from "./modules/prep_genome.nf"
include { BUILD_INDEX } from "./modules/build_index.nf"


//----------------------------------------------------------------------------//
// prep_resources subworkflow
//----------------------------------------------------------------------------//

workflow prep_resources {

    main:

        FETCH_WHITELIST()
        whitelist = FETCH_WHITELIST.out.whitelist

        if (!params.build_STAR_index && params.prebuilt_STAR_index != null) {
            
            star_index = params.prebuilt_STAR_index

        } else {
            
            PREP_GTF()
            PREP_GENOME()
            BUILD_INDEX(PREP_GENOME.out.genome, PREP_GTF.out.gtf)
            star_index = BUILD_INDEX.out.star_index

        }

    emit:

        whitelist = whitelist
        star_index = star_index
    
}


//----------------------------------------------------------------------------//