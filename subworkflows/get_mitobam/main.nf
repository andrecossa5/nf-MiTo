// get_mitobam subworkflows

// Modules
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../common/modules/merge_R1.nf"
include { MERGE_R2 } from "../common/modules/merge_R2.nf"
include { SOLO } from "../common/modules/Solo.nf"
include { FILTER_MAESTER_BAM } from "./modules/filter_bam.nf"
include { MERGE_ATAC_LANES } from "./modules/atac.nf"
include { ATAC_ADD_BARCODE } from "./modules/atac.nf"
include { BWAMEM2_INDEX } from "./modules/atac.nf"
include { BWAMEM2_MT } from "./modules/atac.nf"

 
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

workflow get_atac_bam {

    take:
        ch_fastq_dirs   // [sample, fastq_dir]
        cell_barcodes   // [sample, barcodes]
        genome          // NUMT-masked genome (FASTA)

    main:
        // Get MT-reads from 10x scATAC libraries (barcode in R2): bwa-mem2 on the masked genome
        MERGE_ATAC_LANES(ch_fastq_dirs)
        ATAC_ADD_BARCODE(MERGE_ATAC_LANES.out.reads.combine(cell_barcodes, by:0))
        BWAMEM2_INDEX(genome)
        BWAMEM2_MT(ATAC_ADD_BARCODE.out.reads, BWAMEM2_INDEX.out.index)

    emit:
        mitobam = BWAMEM2_MT.out.bam.combine(cell_barcodes, by:0)
}
//----------------------------------------------------------------------------//
