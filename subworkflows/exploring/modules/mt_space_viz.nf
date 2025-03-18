// VIZ_MT_SPACE module

nextflow.enable.dsl = 2 

//

process VIZ_MT_SPACE {

    tag "${sample}: explore ${job_id}"
    publishDir "${params.output_folder}/${sample}", mode: 'copy'

    input:
    tuple val(job_id), val(sample), val(ch_matrix)
 
    output:
    tuple val(job_id), val(sample), path("${job_id}"), emit: plots

    // Handle CLI from params-file
    def covariate = params.covariate ? "--covariate ${params.covariate}" : ""
    
    script:
    """
    python ${baseDir}/bin/pp/explore_mt_space.py \
    --path_afm ${ch_matrix} \
    --path_tuning ${params.path_tuning} \
    --job_id ${job_id} \
    --sample ${sample} \
    --ncores ${task.cpus} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.path_REDIdb} \
    ${covariate}
    """

    stub:
    """
    mkdir ${job_id}
    cd ${job_id}
    touch aa
    touch bb
    """

}