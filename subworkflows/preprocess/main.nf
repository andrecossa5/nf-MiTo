// Preprocessing subworkflows
nextflow.enable.dsl = 2

// Preprocessing: 10x and MAESTER, raw sequencing data
include { tenx } from "../tenx/main"
include { maester } from "../maester/main" 

// Pre-processing (only MAESTER data), from .bam
include { mitobam } from "../from_bam/main" 

//


// Util to create preprocessing channel
def createPreprocessingChannel() {

    if (params.raw_data_input_type == "fastq") {

        // From fastq, unaligned raw reads
        ch = Channel.fromPath(params.raw_data_input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.fastq_folder, row.library ] }
        
    } else if (params.raw_data_input_type == "bam") {

        // From aligned MT-reads, and a .txt file of valid 10x barcodes
        ch = Channel.fromPath(params.raw_data_input)
            .splitCsv(header: true)
            .map { row -> [ row.sample, row.bam, row.cell_barcodes ] }
    }
    else {
        error "Unsupported raw_data_input_type: ${params.raw_data_input_type}. Available: 'fastq' or 'bam'."
    }

    return ch
}


//


// Raw reads pre-processing
workflow preprocess {

    take:
        ch_preprocessing

    main:
    
        if (params.raw_data_input_type == "fastq") {

            // From raw reads, unaligned
            tenx(ch_preprocessing.filter { it -> it[2] == 'TENX' })
            maester(
                ch_preprocessing.filter { it -> it[2] == 'MAESTER' },
                tenx.out.cell_barcodes,
                tenx.out.bam
            )
            afm = maester.out.afm

        } else if (params.raw_data_input_type == "bam") {

            // From aligned MT-reads, and a .txt file of valid 10x barcodes
            mitobam(ch_preprocessing)
            afm = mitobam.out.afm

        }

    emit:

        afm = Channel.of('MiTo.tree').combine(afm)

}


//