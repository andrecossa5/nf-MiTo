// Preprocessing subworkflows
nextflow.enable.dsl = 2

// Subworkflows
include { tenx } from "../tenx/main"
include { get_tenx_maester_bam } from "../get_mitobam/main" 
include { get_maester_bam } from "../get_mitobam/main" 
include { process_mitobam } from "../process_mitobam/main" 


//


// Util to create preprocessing channel
def createPreprocessingChannel() {

    if (!params.raw_data_input) {
        error "Error: --raw_data_input is required for preprocessing"
    }

    if (params.raw_data_input_type == "fastq") {

        // From raw reads, unaligned
        ch = Channel.fromPath(params.raw_data_input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.fastq_folder, row.library ] }
        
    } else if (params.raw_data_input_type == "fastq, MAESTER") {

        // From raw reads, unaligned (MAESTER) and a .txt file of valid 10x barcodes
        ch = Channel.fromPath(params.raw_data_input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.fastq_folder, row.cell_barcodes ] }
        
    } else if (params.raw_data_input_type == "mitobam") {

        // From aligned MT-reads, and a .txt file of valid 10x barcodes
        ch = Channel.fromPath(params.raw_data_input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.bam, row.cell_barcodes ] }
    }
    else {
        error "Unsupported raw_data_input_type: ${params.raw_data_input_type}. Available: \"fastq\", \"fastq, MAESTER\", or \"mitobam\"."
    }

    return ch
}


//


// Raw data preprocessing, in different flavours
workflow preprocess {

    take:
        ch

    main:
 
        if (params.raw_data_input_type == "fastq") {

            // All 10x and MAESTER MT reads
            tenx_fastqs = ch.filter{it->it[2]=='TENX'}.map{it->tuple(it[0],it[1])}
            maester_fastqs = ch.filter{it->it[2]=='MAESTER'}.map{it->tuple(it[0],it[1])}
            tenx(tenx_fastqs)
            get_tenx_maester_bam(maester_fastqs, tenx.out.cell_barcodes_QC, tenx.out.bam)
            ch_mitobam = get_tenx_maester_bam.out.mitobam

        } else if (params.raw_data_input_type == "fastq, MAESTER") {

            // MAESTER MT reads
            maester_fastqs = ch.map{it->tuple(it[0],it[1])}
            cell_barcodes = ch.map{it->tuple(it[0],it[2])}
            get_maester_bam(maester_fastqs, cell_barcodes)
            ch_mitobam = get_maester_bam.out.mitobam

        } else if (params.raw_data_input_type == "mitobam") {

            // Previously aligned (and filtered for chrM) MT-reads
            ch_mitobam = ch

        }
     
        // Process mitobam, for a list of target cell barcodes
        process_mitobam(ch_mitobam)
        afm = process_mitobam.out.afm
        
    emit:

        afm = Channel.of('MiTo.tree').combine(afm)

}


//