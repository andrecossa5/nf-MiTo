// mitobam subworkflow

// Include here
nextflow.enable.dsl = 2
include { SPLIT_BARCODES } from "./modules/split_barcodes.nf"
include { FILTER_BAM_CB } from "./modules/filter_bam_cb.nf"
include { SPLIT_BAM } from "./modules/split_bam.nf"
include { EXTRACT_FASTA } from "./modules/extract_fasta.nf"
include { MAEGATK } from "./modules/maegatk.nf"
include { COLLAPSE_MAEGATK } from "./modules/maegatk.nf"

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
// mitobam subworkflow
//----------------------------------------------------------------------------//

workflow process_mitobam {
     
    take:
        ch_mitobam
        genome

    main:

        EXTRACT_FASTA(params.string_MT, genome)
        SPLIT_BARCODES(ch_mitobam)
        ch_chunks = SPLIT_BARCODES.out.bam_chunks.flatMap { 
            sample_name, bam, file_paths ->
            if (file_paths instanceof List) {
                file_paths.collect { file_path -> [sample_name, bam, file_path] }
            } else {
                [[sample_name, bam, file_paths]]
            }
        }
        FILTER_BAM_CB(ch_chunks)

        if (params.pp_method == "maegatk") {

            SPLIT_BAM(FILTER_BAM_CB.out.bam)
            ch_cell_bams = processCellBams(SPLIT_BAM.out.cell_bams)
            MAEGATK(ch_cell_bams, EXTRACT_FASTA.out.fasta)
            COLLAPSE_MAEGATK(MAEGATK.out.tables.groupTuple(by:0))
            ch_output = COLLAPSE_MAEGATK.out.ch_output

        } else {

            error "Unsupported pp_method: ${params.pp_method}. Available here: maegatk."

        }

    emit:

        tables = ch_output

}

//----------------------------------------------------------------------------//