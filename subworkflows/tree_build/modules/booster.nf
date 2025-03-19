// BOOSTER module

nextflow.enable.dsl = 2


process BOOSTER {

    tag "${sample}"

    input:
    tuple val(job_id),
        val(sample), 
        path(trees)

    output:
    tuple val(job_id),
        val(sample), 
        path("final_tree.newick"), emit: final_tree
     
    script:
    """
    mv rep_observed.newick observed.newick
    booster -i observed.newick -b rep* -o final_tree.newick -@ ${task.cpus} -a ${params.support_method}
    """

    stub:
    """
    touch final_tree.newick
    """

}
