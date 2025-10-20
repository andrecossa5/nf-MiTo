// ONESAMPLE module

nextflow.enable.dsl = 2 

//

process ONESAMPLE {

    tag "${sample}: tuning${job_id}"
    label 'MiTo'

    input:
    tuple val(sample), 
        path(path_afm), 
        val(min_n_positive),
        val(af_confident_detection),
        val(min_n_confidently_detected),
        val(min_mean_AD_in_positives),
        val(min_AD),
        val(bin_method),
        val(t_prob),
        val(min_cell_prevalence),
        val(job_id)
 
    output:
    tuple val(job_id), 
        val(sample), 
        path("job_${job_id}_metrics.csv"),
        path("job_${job_id}_options.csv"), emit: stats
    
    script:
    
    // Handle CLI from params (only conditional for null-defaulting parameters)
    def lineage_column = params.lineage_column ? "--lineage_column ${params.lineage_column}" : ""

    """
    python ${baseDir}/bin/afm_preprocess/onesample.py \
    --path_afm ${path_afm} \
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
    --min_n_positive ${min_n_positive} \
    --af_confident_detection ${af_confident_detection} \
    --min_n_confidently_detected ${min_n_confidently_detected} \
    --min_mean_AD_in_positives ${min_mean_AD_in_positives} \
    --min_mean_DP_in_positives ${params.min_mean_DP_in_positives} \
    --t_prob ${t_prob} \
    --t_vanilla ${params.t_vanilla} \
    --min_AD ${min_AD} \
    --min_cell_prevalence ${min_cell_prevalence} \
    --bin_method ${bin_method} \
    --k ${params.k} \
    --gamma ${params.gamma} \
    --min_n_var ${params.min_n_var} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.distance_metric} \
    --filter_dbs ${params.filter_dbs} \
    --filter_moran ${params.filter_moran} \
    --ncores ${task.cpus} \
    --max_fraction_unassigned ${params.max_fraction_unassigned} \
    ${lineage_column}
    """

    stub:
    """
    touch job_${job_id}_metrics.csv
    touch job_${job_id}_options.csv
    """

}