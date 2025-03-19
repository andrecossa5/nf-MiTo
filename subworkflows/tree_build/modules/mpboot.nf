// MPBOOT module

nextflow.enable.dsl = 2

process MPBOOT {

    tag "${sample}: ${job_id}, rep=${rep}"

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
    python ${baseDir}/bin/build_tree/create_fasta.py ${afm}
    mpboot -s genotypes.fa
    mv genotypes.fa.treefile rep_${rep}.newick
    """

    stub:
    """
    touch rep_${rep}.newick
    """

}
