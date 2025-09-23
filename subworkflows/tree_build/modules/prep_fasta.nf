// PREP_FASTA module

nextflow.enable.dsl = 2

process PREP_FASTA {

    tag "${sample}: ${job_id}, rep=${rep}"
    label 'MiTo'

    input:
    tuple val(job_id),
        val(sample), 
        val(rep),
        val(afm)

    output:
    tuple val(job_id),
        val(sample), 
        val(rep),
        path("genotypes.fa"), emit: fasta
    
    script:
    """
    python ${baseDir}/bin/tree_build/create_fasta.py ${afm}
    """

    stub:
    """
    touch genotypes.fa
    """

}