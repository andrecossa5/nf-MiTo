// ONESAMPLE module

nextflow.enable.dsl = 2 

//

process SUMMARISE {

    publishDir "${params.output_folder}", mode: 'copy'

    input:
    path(files)
 
    output:
    tuple path("all_metrics_final.csv"), path("all_options_final.csv"), emit: summary

    script: 
    """
    (echo "sample,job_id,option,value"; cat *_options.csv) > all_options_final.csv
    (echo "sample,job_id,metric,value"; cat *_metrics.csv) > all_metrics_final.csv
    """

    stub:
    """
    touch all_options_final.csv
    touch all_metrics_final.csv
    """

}