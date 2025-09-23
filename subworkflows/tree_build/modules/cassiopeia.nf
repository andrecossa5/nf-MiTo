// CASSIOPEIA module

nextflow.enable.dsl = 2


process CASSIOPEIA {

    tag "${sample}: ${job_id}, rep=${rep}"
    label 'MiTo'

    input:
    tuple val(job_id),
        val(sample), 
        val(rep),
        path(afm)

    output:
    tuple val(job_id),
        val(sample), 
        path("*.newick"), emit: tree
    
    script:
    // Handle CLI from params (only conditional for null-defaulting parameters)
    def path_tuning = params.path_tuning ? "--path_tuning ${params.path_tuning}" : ""

    """
    python ${baseDir}/bin/tree_build/build_cassiopeia.py \
    --path_afm ${afm} \
    --sample ${sample} \
    --job_id ${job_id} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.distance_metric} \
    --boot_replicate ${rep} \
    --ncores ${task.cpus} \
    ${path_tuning}
    """

    stub:
    """
    touch rep_${rep}.newick
    """

}


//