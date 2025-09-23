// IQTREE module

nextflow.enable.dsl = 2

process IQTREE {

    tag "${sample}: ${job_id}, rep=${rep}"
    label 'phylo'

    input:
    tuple val(job_id),
        val(sample), 
        val(rep),
        path(genotypes)

    output:
    tuple val(job_id),
        val(sample), 
        path("*.newick"), emit: tree
    
    script:
    """
    iqtree -s ${genotypes} -m GTR
    mv ${genotypes}.treefile rep_${rep}.newick
    """

    stub:
    """
    touch rep_${rep}.newick
    """

}