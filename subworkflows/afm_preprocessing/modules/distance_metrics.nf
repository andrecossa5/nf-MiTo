// DISTANCES module

nextflow.enable.dsl = 2 

//

process DISTANCE_METRICS {

    tag "${sample}: ${job_id}"
    publishDir "${params.output_folder}/${sample}/${job_id}", mode: 'copy'

    input:
    tuple val(job_id),
        val(sample), 
        val(rep),
        val(afm)

    output:
    tuple val(job_id),
        val(sample),  
        path("afm_filtered.h5ad"), emit: distance_metrics
    
    // Handle CLI args
    def lineage_column = params.lineage_column ? "--lineage_column ${params.lineage_column}" : ""
    def K = params.K ? "--K ${params.K}" : ""

    script:
    """
    python ${baseDir}/bin/pp/distance_metrics.py \
    --afm "${afm}" \
    --replicates "${rep}" \
    --job_id ${job_id} \
    ${K} \
    ${lineage_column}
    """

    stub:
    """
    touch afm_filtered.h5ad
    """

}

