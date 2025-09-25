// prep_resources workflow

// Include here
nextflow.enable.dsl = 2
include { FETCH_WHITELIST } from "./modules/fetch_whitelist.nf"
include { PREP_GTF } from "./modules/prep_gtf.nf"
include { PREP_GENOME } from "./modules/prep_genome.nf"
include { BUILD_INDEX } from "./modules/build_index.nf"
include { VALIDATE_STAR_INDEX } from "./modules/validate_star_index.nf"


//----------------------------------------------------------------------------//
// prep_resources and prep_genome subworkflows
//----------------------------------------------------------------------------//

workflow prep_resources {

    main:

        FETCH_WHITELIST()
        whitelist = FETCH_WHITELIST.out.whitelist

        if (!params.build_STAR_index && params.prebuilt_STAR_index != null) {
            
            VALIDATE_STAR_INDEX(params.prebuilt_STAR_index)
            star_index = VALIDATE_STAR_INDEX.out.star_index
            genome = VALIDATE_STAR_INDEX.out.genome

        } else {
            
            PREP_GTF()
            PREP_GENOME()
            BUILD_INDEX(PREP_GENOME.out.genome, PREP_GTF.out.gtf)
            star_index = BUILD_INDEX.out.star_index
            genome = PREP_GENOME.out.genome

        }

    emit:

        whitelist = whitelist
        star_index = star_index
        genome = genome
    
}


//


workflow prep_genome {

    main:

        PREP_GENOME()

    emit:

        genome = PREP_GENOME.out.genome
    
}


//


//----------------------------------------------------------------------------//