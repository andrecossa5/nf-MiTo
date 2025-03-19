 // TREE_METRICS module

nextflow.enable.dsl = 2


process TREE_METRICS {

    tag "${sample}: ${job_id}"
    publishDir "${params.output_folder}/${sample}/${job_id}", mode: 'copy'

    input:
    tuple val(job_id),
        val(sample), 
        path(tree)

    output:
    tuple val(job_id),
        val(sample), 
        path("tree_metrics.csv"), emit: metrics

    // Handle CLI args
    def lineage_column = params.lineage_column ? "--lineage_column ${params.lineage_column}" : ""
    
    script:
    """
    python ${baseDir}/bin/tree_process/tree_metrics.py \
    --tree ${tree} \
    ${lineage_column} \
    --job_id ${job_id}
    """

    stub:
    """
    touch tree_metrics.csv
    """

}
