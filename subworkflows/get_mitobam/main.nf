// get_mitobam subworkflows

// Modules
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../common/modules/merge_R1.nf"
include { MERGE_R2 } from "../common/modules/merge_R2.nf"
include { SOLO } from "../common/modules/Solo.nf"
include { FILTER_10X_BAM } from "./modules/filter_bam.nf"
include { FILTER_MAESTER_BAM } from "./modules/filter_bam.nf"
include { MERGE_BAM } from "./modules/merge_bams.nf"

 
// 

def processCellBams(cell_bams) {
    return cell_bams
        .map { it ->
            def sample = it[0]
            def paths = it[1]
            return paths.collect { cell_path ->
                def path_splitted = cell_path.toString().split('/')
                def cell = path_splitted[-1].toString().split('\\.')[0]
                return [sample, cell, cell_path]
            }
        }
        .flatMap { it }
}

//

//----------------------------------------------------------------------------//
// get_tenx_maester_bam, get_maester_bam subworkflows
//----------------------------------------------------------------------------//

workflow get_tenx_maester_bam {
     
    take:
        ch_input
        cell_barcodes
        tenx_bam  

    main:

        // Get MT-reads from 10x and MAESTER libraries
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))
        FILTER_MAESTER_BAM(SOLO.out.bam)
        FILTER_10X_BAM(tenx_bam)
        MERGE_BAM(FILTER_10X_BAM.out.bam.combine(FILTER_MAESTER_BAM.out.bam, by:0))

    emit:

        mitobam = MERGE_BAM.out.bam.combine(cell_barcodes, by:0)

}


//


workflow get_maester_bam {
     
    take:
        ch_input
        cell_barcodes

    main:

        // Get MT-reads from MAESTER libraries
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))
        FILTER_MAESTER_BAM(SOLO.out.bam)

    emit:

        mitobam = FILTER_MAESTER_BAM.out.bam.combine(cell_barcodes, by:0)

}

//----------------------------------------------------------------------------//