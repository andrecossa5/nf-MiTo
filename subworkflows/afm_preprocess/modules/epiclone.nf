// EPICLONE module

nextflow.enable.dsl = 2 

//

process EPICLONE {

    tag "${sample}: ${job_id}"
    label 'MiTo'

    input:
    tuple val(job_id), val(sample), val(ch_matrix)
 
    output:
    tuple val(job_id), val(sample), path("afm.h5ad"), emit: afm
    
    script:
    """
    python ${baseDir}/bin/afm_preprocess/epiclone.py \
    --path_afm ${ch_matrix} \
    --metric ${params.distance_metric} \
    --n_cores ${task.cpus}
    """

    stub:
    """
    touch afm.h5ad
    """

}
