// get_mitobam subworkflows

// Modules
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../common/modules/merge_R1.nf"
include { MERGE_R2 } from "../common/modules/merge_R2.nf"
include { SOLO } from "../common/modules/Solo.nf"
include { FILTER_MAESTER_BAM } from "./modules/filter_bam.nf"

 
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
// get_maester_bam subworkflow
//----------------------------------------------------------------------------//

workflow get_maester_bam {
     
    take:
        ch_fastqs
        cell_barcodes
        whitelist
        star_index

    main:

        // Get MT-reads from MAESTER libraries
        MERGE_R1(ch_fastqs)
        MERGE_R2(ch_fastqs)
        fastqs = MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0)
        SOLO(fastqs, whitelist, star_index)
        FILTER_MAESTER_BAM(SOLO.out.bam)

    emit:

        mitobam = FILTER_MAESTER_BAM.out.bam.combine(cell_barcodes, by:0)

}

//----------------------------------------------------------------------------//