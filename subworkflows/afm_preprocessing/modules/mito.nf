// MITO module

nextflow.enable.dsl = 2 

//

process MITO {

    tag "${sample}: ${job_id}"

    input:
    tuple val(job_id), val(sample), val(ch_matrix)
 
    output:
    tuple val(job_id), val(sample), path("afm.h5ad"), emit: afm

    // Handle CLI from params-file
    def path_tuning = params.path_tuning ? "--path_tuning ${params.path_tuning}" : ""
    def cell_filter = params.cell_filter ? "--cell_filter ${params.cell_filter}" : ""
    def filtering = params.filtering ? "--filtering ${params.filtering}" : ""
    def min_cell_number = params.min_cell_number ? "--min_cell_number ${params.min_cell_number}" : ""
    def min_cov = params.min_cov ? "--min_cov ${params.min_cov}" : ""
    def min_var_quality = params.min_var_quality ? "--min_var_quality ${params.min_var_quality}" : ""
    def min_frac_negative = params.min_frac_negative ? "--min_frac_negative ${params.min_frac_negative}" : ""
    def min_n_positive = params.min_n_positive ? "--min_n_positive ${params.min_n_positive}" : ""
    def af_confident_detection = params.af_confident_detection ? "--af_confident_detection ${params.af_confident_detection}" : ""
    def min_n_confidently_detected = params.min_n_confidently_detected ? "--min_n_confidently_detected ${params.min_n_confidently_detected}" : ""
    def min_mean_AD_in_positives = params.min_mean_AD_in_positives ? "--min_mean_AD_in_positives ${params.min_mean_AD_in_positives}" : "" 
    def min_mean_DP_in_positives = params.min_mean_DP_in_positives ? "--min_mean_DP_in_positives ${params.min_mean_DP_in_positives}" : ""
    def lineage_column = params.lineage_column ? "--lineage_column ${params.lineage_column}" : ""
    def t_prob = params.t_prob ? "--t_prob ${params.t_prob}" : ""
    def min_AD = params.min_AD ? "--min_AD ${params.min_AD}" : ""
    def t_vanilla = params.t_vanilla ? "--t_vanilla ${params.t_vanilla}" : ""
    def min_cell_prevalence = params.min_cell_prevalence ? "--min_cell_prevalence ${params.min_cell_prevalence}" : ""
    def bin_method = params.bin_method ? "--bin_method ${params.bin_method}" : ""
    def cassiopeia_solver = params.cassiopeia_solver ? "--solver ${params.cassiopeia_solver}" : ""
    def distance_metric = params.distance_metric ? "--metric ${params.distance_metric}" : ""
    def k = params.k ? "--k ${params.k}" : ""
    def gamma = params.gamma ? "--gamma ${params.gamma}" : ""
    def min_n_var = params.min_n_var ? "--min_n_var ${params.min_n_var}" : ""
    
    
    script:
    """
    python ${baseDir}/bin/pp/MiTo.py \
    --path_afm ${ch_matrix} \
    ${path_tuning} \
    --job_id ${job_id} \
    --sample ${sample} \
    ${cell_filter} \
    ${filtering} \
    ${min_cell_number} \
    ${min_cov} \
    ${min_var_quality} \
    ${min_frac_negative} \
    ${min_n_positive} \
    ${af_confident_detection} \
    ${min_n_confidently_detected} \
    ${min_mean_AD_in_positives} \
    ${min_mean_DP_in_positives} \
    ${t_prob} \
    ${t_vanilla} \
    ${min_AD} \
    ${min_cell_prevalence} \
    ${bin_method} \
    ${lineage_column} \
    ${cassiopeia_solver} \
    ${distance_metric} \
    ${k} \
    ${gamma} \
    ${min_n_var} \
    --ncores ${task.cpus} \
    --path_dbSNP ${params.path_dbSNP} \
    --path_REDIdb ${params.path_REDIdb}
    """

    stub:
    """
    touch afm.h5ad
    """

}
