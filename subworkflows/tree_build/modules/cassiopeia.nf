// CASSIOPEIA module

nextflow.enable.dsl = 2


process CASSIOPEIA {

    tag "${sample}: ${job_id}, rep=${rep}"

    input:
    tuple val(job_id),
        val(sample), 
        val(rep),
        path(afm)

    output:
    tuple val(job_id),
        val(sample), 
        path("*.newick"), emit: tree
    
    // Handle CLI args
    def path_tuning = params.path_tuning ? "--path_tuning ${params.path_tuning}" : ""

    script:
    """
    python ${baseDir}/bin/tree_build/build_cassiopeia.py \
    --path_afm ${afm} \
    ${path_tuning} \
    --sample ${sample} \
    --job_id ${job_id} \
    --solver ${params.cassiopeia_solver} \
    --metric ${params.distance_metric} \
    --boot_replicate ${rep} \
    --ncores ${task.cpus}
    """

    stub:
    """
    touch rep_${rep}.newick
    """

}


//