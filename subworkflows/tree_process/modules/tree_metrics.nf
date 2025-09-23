 // TREE_METRICS module

nextflow.enable.dsl = 2


process TREE_METRICS {

    tag "${sample}: ${job_id}"
    label 'MiTo'
    publishDir "${params.output_folder}/${sample}/${job_id}", mode: 'copy'

    input:
    tuple val(job_id),
        val(sample), 
        path(tree)

    output:
    tuple val(job_id),
        val(sample), 
        path("tree_metrics.csv"), emit: metrics

    script:
    // Handle CLI from params (only conditional for null-defaulting parameters)
    def lineage_column = params.lineage_column ? "--lineage_column ${params.lineage_column}" : ""
    
    """
    python ${baseDir}/bin/tree_process/tree_metrics.py \
    --tree ${tree} \
    --job_id ${job_id} \
    ${lineage_column}
    """

    stub:
    """
    touch tree_metrics.csv
    """

}
