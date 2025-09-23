// PATH module

nextflow.enable.dsl = 2


process PATH {

    tag "${sample}: ${job_id}"
    label 'phylo'
    publishDir "${params.output_folder}/${sample}/${job_id}", mode: 'copy'

    input:
    tuple val(job_id),
        val(sample), 
        path(tree)

    output:
    tuple val(job_id),
        val(sample), 
        path("phylocorr.csv"), emit: phylocorr
    
    script:
    """
    Rscript ${baseDir}/bin/tree_process/PATH.r ${tree} ${params.path_meta} ${params.lineage_column}
    """

    stub:
    """
    touch phylocorr.csv
    """

}
