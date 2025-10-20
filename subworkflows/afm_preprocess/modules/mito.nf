// MITO module

nextflow.enable.dsl = 2 

//

process MITO {

    tag "${sample}: ${job_id}"
    label 'MiTo'

    input:
    tuple val(job_id), val(sample), val(ch_matrix)
 
    output:
    tuple val(job_id), val(sample), path("afm.h5ad"), emit: afm

    // Handle CLI from params-file (only conditional for null-defaulting parameters)
    def path_tuning = params.path_tuning ? "--path_tuning ${params.path_tuning}" : ""
    def lineage_column = params.lineage_column ? "--lineage_column ${params.lineage_column}" : ""
    
    script:
    """
    python ${baseDir}/bin/afm_preprocess/MiTo.py \
    --path_afm ${ch_matrix} \
    --job_id ${job_id} \
    --sample ${sample} \
    --cell_filter ${params.cell_filter} \
    --mean_cov_all ${params.mean_cov_all} \
    --median_cov_target ${params.median_cov_target} \
    --min_perc_covered_sites ${params.min_perc_covered_sites} \
    --filtering ${params.filtering} \
    --min_cell_number ${params.min_cell_number} \
    --min_cov ${params.min_cov} \
    --min_var_quality ${params.min_var_quality} \
    --min_frac_negative ${params.min_frac_negative} \
    --min_n_positive ${params.min_n_positive} \
    --af_confident_detection ${params.af_confident_detection} \
    --min_n_confidently_detected ${params.min_n_confidently_detected} \
    --min_mean_AD_in_positives ${params.min_mean_AD_in_positives} \
    --min_mean_DP_in_positives ${params.min_mean_DP_in_positives} \
    --t_prob ${params.t_prob} \
    --t_vanilla ${params.t_vanilla} \
    --min_AD ${params.min_AD} \
    --min_cell_prevalence ${params.min_cell_prevalence} \
    --bin_method ${params.bin_method} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.distance_metric} \
    --k ${params.k} \
    --gamma ${params.gamma} \
    --min_n_var ${params.min_n_var} \
    --filter_moran ${params.filter_moran} \
    --filter_dbs ${params.filter_dbs} \
    --spatial_metrics ${params.spatial_metrics} \
    --ncores ${task.cpus} \
    ${path_tuning} \
    ${lineage_column} \
    """

    stub:
    """
    touch afm.h5ad
    """

}
