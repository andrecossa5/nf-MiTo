// FILTER_MITOBAM module

nextflow.enable.dsl = 2

//

process FILTER_BAM_CB {

    tag "${sample_name}"
    label 'scLT'
    
    input:
    tuple val(sample_name), path(bam), path(CBs)
    
    output:
    tuple val(sample_name), path("*.bam"), emit: bam

    script:
    """
    python ${baseDir}/bin/preprocess/filter_bam_cb.py ${bam} ${CBs}
    """

    stub:
    """
    id=\$(basename ${CBs} | cut -d '_' -f2 | cut -d '.' -f1)
    touch "filtered_\${id}.bam"
    """

} 
