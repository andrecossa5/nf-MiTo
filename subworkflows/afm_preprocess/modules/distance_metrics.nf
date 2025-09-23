// DISTANCES module

nextflow.enable.dsl = 2 

//

process DISTANCE_METRICS {

    tag "${sample}: ${job_id}"
    label 'MiTo'
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
    
    script:

    // Handle CLI from params (only conditional for null-defaulting parameters)
    def lineage_column = params.lineage_column ? "--lineage_column ${params.lineage_column}" : ""

    """
    python ${baseDir}/bin/afm_preprocess/distance_metrics.py \
    --afm "${afm}" \
    --replicates "${rep}" \
    --job_id ${job_id} \
    --K ${params.K} \
    ${lineage_column}
    """

    stub:
    """
    touch afm_filtered.h5ad
    """

}

