// PROCESS_TREE module

nextflow.enable.dsl = 2


process PROCESS_TREE {

    tag "${sample}: ${job_id}"
    label 'MiTo'
    publishDir "${params.output_folder}/${sample}/${job_id}", mode: 'copy'

    input:
    tuple val(job_id),
        val(sample), 
        val(observed),
        path(afm),
        path(tree)

    output:
    tuple val(job_id),
        val(sample), 
        path("annotated_tree.pickle"), emit: annotated_tree

    script:
    """
    python ${baseDir}/bin/tree_process/annotate_tree.py \
    --tree ${tree} \
    --afm ${afm} \
    --annotate_tree ${params.annotate_tree} \
    --max_fraction_unassigned ${params.max_fraction_unassigned}
    """

    stub:
    """
    touch annotated_tree.pickle
    """

}
