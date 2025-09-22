// VIZ_MT_SPACE module

nextflow.enable.dsl = 2 

//

process VIZ_MT_SPACE {

    tag "${sample}: explore ${job_id}"
    label 'MiTo'
    publishDir "${params.output_folder}/${sample}", mode: 'copy'

    input:
    tuple val(job_id), val(sample), val(ch_matrix)
 
    output:
    tuple val(job_id), val(sample), path("${job_id}"), emit: plots

    // Handle CLI from params-file
    def covariate = params.covariate ? "--covariate ${params.covariate}" : ""
    def coverage_input = params.coverage_input ? params.coverage_input : ""
    
    script:
    """
    python ${baseDir}/bin/afm_preprocess/explore_mt_space.py \
    --path_afm ${ch_matrix} \
    --path_tuning ${params.path_tuning} \
    --job_id ${job_id} \
    --sample ${sample} \
    --coverage_input ${coverage_input} \
    --ncores ${task.cpus} \
    --filter_dbs ${params.filter_dbs} \
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