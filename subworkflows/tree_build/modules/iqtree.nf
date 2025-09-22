// IQTREE module

nextflow.enable.dsl = 2

process IQTREE {

    tag "${sample}: ${job_id}, rep=${rep}"
    label 'phylo'

    input:
    tuple val(job_id),
        val(sample), 
        val(rep),
        val(afm)

    output:
    tuple val(job_id),
        val(sample), 
        path("*.newick"), emit: tree
    
    script:
    """
    python ${baseDir}/bin/tree_build/create_fasta.py ${afm}
    iqtree -s genotypes.fa -m GTR
    mv genotypes.fa.treefile rep_${rep}.newick
    """

    stub:
    """
    touch rep_${rep}.newick
    """

}