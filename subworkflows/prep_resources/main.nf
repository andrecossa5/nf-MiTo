// prep_resources workflow

// Include here
nextflow.enable.dsl = 2
include { FETCH_WHITELIST } from "./modules/fetch_whitelist.nf"
include { PREP_GTF } from "./modules/prep_gtf.nf"
include { BUILD_INDEX } from "./modules/build_index.nf"
include { VALIDATE_INDEX } from "./modules/validate_index.nf"


// The reference is given by the user: --masked_ref, a NUMT-masked, uncompressed genome FASTA
def masked_ref() {
    return Channel.value(file(params.masked_ref, checkIfExists: true))
}

//----------------------------------------------------------------------------//
// prep_resources and prep_genome subworkflows
//----------------------------------------------------------------------------//

workflow prep_resources {

    main:

        FETCH_WHITELIST()
        whitelist = FETCH_WHITELIST.out.whitelist

        if (!params.build_STAR_index && params.prebuilt_STAR_index != null) {
            
            VALIDATE_INDEX()
            star_index = VALIDATE_INDEX.out.star_index

        } else {
            
            PREP_GTF()
            BUILD_INDEX(masked_ref(), PREP_GTF.out.gtf)
            star_index = BUILD_INDEX.out.star_index

        }

    emit:

        whitelist = whitelist
        star_index = star_index
        genome = masked_ref()
    
}


//


workflow prep_genome {

    main:

        genome = masked_ref()

    emit:

        genome = genome
    
}


//


//----------------------------------------------------------------------------//