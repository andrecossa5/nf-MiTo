// LEIDEN module

nextflow.enable.dsl = 2 

//

process LEIDEN {

    tag "${sample}: ${job_id}"
    publishDir "${params.output_folder}/${sample}", mode: 'copy'
    label 'MiTo'

    input:
    tuple val(job_id), val(sample), val(ch_matrix)
 
    output:
    tuple val(job_id), val(sample), path("${job_id}_leiden.pickle"), emit: results
    
    script:
    """
    python ${baseDir}/bin/_bench/leiden.py \
    --path_afm ${ch_matrix} \
    --sample ${sample} \
    --job_id ${job_id}
    """

    stub:
    """
    touch "${job_id}_leiden.pickle"
    """

}
