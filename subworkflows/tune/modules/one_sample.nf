// ONESAMPLE module

nextflow.enable.dsl = 2 

//

process ONESAMPLE {

    tag "${sample}: tuning${job_id}"

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
    
    // Hadle CLI
    def cell_filter = params.cell_filter ? "--cell_filter ${params.cell_filter}" : ""
    def filtering = params.filtering ? "--filtering ${params.filtering}" : ""
    def min_cell_number = params.min_cell_number ? "--min_cell_number ${params.min_cell_number}" : ""
    def min_cov = params.min_cov ? "--min_cov ${params.min_cov}" : ""
    def min_var_quality = params.min_var_quality ? "--min_var_quality ${params.min_var_quality}" : ""
    def min_frac_negative = params.min_frac_negative ? "--min_frac_negative ${params.min_frac_negative}" : ""
    def min_mean_DP_in_positives = params.min_mean_DP_in_positives ? "--min_mean_DP_in_positives ${params.min_mean_DP_in_positives}" : ""
    def lineage_column = params.lineage_column ? "--lineage_column ${params.lineage_column}" : ""
    def k = params.k ? "--k ${params.k}" : ""
    def gamma = params.gamma ? "--gamma ${params.gamma}" : ""
    def min_n_var = params.min_n_var ? "--min_n_var ${params.min_n_var}" : ""
    def max_fraction_unassigned = params.max_fraction_unassigned ? "--max_fraction_unassigned ${params.max_fraction_unassigned}" : ""

    script: 
    """
    python ${baseDir}/bin/afm_preprocess/onesample.py \
    --path_afm ${path_afm} \
    --job_id ${job_id} \
    --sample ${sample} \
    ${cell_filter} \
    ${filtering} \
    ${min_cell_number} \
    ${min_cov} \
    ${min_var_quality} \
    ${min_frac_negative} \
    --min_n_positive ${min_n_positive} \
    --af_confident_detection ${af_confident_detection} \
    --min_n_confidently_detected ${min_n_confidently_detected} \
    --min_mean_AD_in_positives ${min_mean_AD_in_positives} \
    ${min_mean_DP_in_positives} \
    --t_prob ${t_prob} \
    --t_vanilla ${params.t_vanilla} \
    --min_AD ${min_AD} \
    --min_cell_prevalence ${min_cell_prevalence} \
    --bin_method ${bin_method} \
    ${k} \
    ${gamma} \
    ${min_n_var} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.distance_metric} \
    ${lineage_column} \
    --ncores ${task.cpus} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.path_REDIdb} \
    ${max_fraction_unassigned}
    """

    stub:
    """
    touch job_${job_id}_metrics.csv
    touch job_${job_id}_options.csv
    """

}