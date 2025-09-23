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

    script:
    
    // Handle CLI from params (only conditional for null-defaulting parameters)
    def covariate = params.covariate ? "--covariate ${params.covariate}" : ""
    def coverage_input = params.coverage_input ? "--coverage_input ${params.coverage_input}" : ""
    def path_tuning = params.path_tuning ? "--path_tuning ${params.path_tuning}" : ""

    """
    python ${baseDir}/bin/afm_preprocess/explore_mt_space.py \
    --path_afm ${ch_matrix} \
    --job_id ${job_id} \
    --sample ${sample} \
    --ncores ${task.cpus} \
    --filter_dbs ${params.filter_dbs} \
    --max_fraction_unassigned ${params.max_fraction_unassigned} \
    ${path_tuning} \
    ${coverage_input} \
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