// sc_pp workflow

// Include here
nextflow.enable.dsl = 2
include { MERGE_R1 } from "../common/modules/merge_R1.nf"
include { MERGE_R2 } from "../common/modules/merge_R2.nf"
include { SOLO } from "../common/modules/Solo.nf"

// 

process publish_tenx {

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
    mv raw Solo
    mv filtered Solo
    mv Features.stats Solo
    mv Summary.csv Solo
    mv Aligned.sortedByCoord.out.bam Solo
    """

    stub:
    """
    mkdir Solo.output
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
        publish_input = SOLO.out.raw
            .combine(SOLO.out.filtered, by:0)
            .combine(SOLO.out.stats, by:0)
            .combine(SOLO.out.summary, by:0)
            .combine(SOLO.out.bam, by:0)
        publish_tenx(publish_input)

    emit:
    
        cell_barcodes = SOLO.out.cell_barcodes
        bam = SOLO.out.bam

}


//----------------------------------------------------------------------------//