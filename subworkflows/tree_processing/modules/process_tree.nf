// PROCESS_TREE module

nextflow.enable.dsl = 2


process PROCESS_TREE {

    tag "${sample}: ${job_id}"
    publishDir "${params.outdir}/${sample}/${job_id}", mode: 'copy'

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

    // Handle CLI args
    def annotate_tree = params.annotate_tree ? "--annotate_tree ${params.annotate_tree}" : ""
    
    script:
    """
    python ${baseDir}/bin/process_tree/annotate_tree.py \
    --tree ${tree} \
    --afm ${afm} \
    ${annotate_tree}
    """

    stub:
    """
    touch annotated_tree.pickle
    """

}
