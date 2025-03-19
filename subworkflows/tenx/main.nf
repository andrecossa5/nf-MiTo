// sc_pp workflow

// Include here
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../common/modules/merge_R1.nf"
include { MERGE_R2 } from "../common/modules/merge_R2.nf"
include { SOLO } from "../common/modules/Solo.nf"
include { QC } from "./modules/QC.nf"

// 

process publish_solo {

    tag "${sample_name}"

    // Publish
    publishDir "${params.output_folder}/${sample_name}/", mode: 'copy'

    input:
    tuple val(sample_name),
          path (raw),
          path (filtered),
          path (stats), 
          path (summary),
          path (bam)

    output:
    path('Solo.output'), emit: solo 

    script:
    """
    mkdir Solo.output
    mv raw Solo.output
    mv filtered Solo.output
    mv Features.stats Solo.output
    mv Summary.csv Solo.output
    mv Aligned.sortedByCoord.out.bam Solo.output
    """

    stub:
    """
    mkdir Solo.output
    cd Solo.output
    touch Features.stats
    touch Summary.csv
    touch Aligned.sortedByCoord.out.bam
    mkdir raw 
    mkdir filtered
    cd filtered 
    touch matrix.mtx.gz
    touch barcodes.txt.gz
    touch features.txt.gz
    """

}

// 


//----------------------------------------------------------------------------//
// tenx subworkflow
//----------------------------------------------------------------------------//


workflow tenx {

    take:
        ch_input

    main:
    
        MERGE_R1(ch_input)
        MERGE_R2(ch_input)
        SOLO(MERGE_R1.out.R1.combine(MERGE_R2.out.R2, by:0))
        ch_solo = SOLO.out.raw
            .combine(SOLO.out.filtered, by:0)
            .combine(SOLO.out.stats, by:0)
            .combine(SOLO.out.summary, by:0)
            .combine(SOLO.out.bam, by:0)
        publish_solo(ch_solo)
        QC(SOLO.out.filtered)

    emit:
    
        gene_expression_matrix = QC.out.gene_expression_matrix
        cell_barcodes_QC = QC.out.cell_barcodes
        bam = SOLO.out.bam

}


//----------------------------------------------------------------------------//