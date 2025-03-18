// DISTANCES module

nextflow.enable.dsl = 2 

//

process DISTANCES {

    tag "${sample}: ${job_id}, rep ${rep}"

    input:
    tuple val(job_id),
        val(sample),  
        path(afm), 
        val(rep)

    output:
    tuple val(job_id),
        val(sample), 
        val(rep), 
        path('afm_new.h5ad'), emit: distances
    
    script:
    """
    python ${baseDir}/bin/pp/distances.py \
    --afm ${afm} \
    --n_cores ${task.cpus} \
    --boot_strategy ${params.boot_strategy} \
    --boot_replicate ${rep} \
    --frac_char_resampling ${params.frac_char_resampling}
    """

    stub:
    """
    touch afm_new.h5ad
    """

}


