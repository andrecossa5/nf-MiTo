// PREP_GTF module

nextflow.enable.dsl = 2

//

process PREP_GTF {

  tag "${sample_name}"
  label 'scLT'

  output:
  path('filtered_genes.gtf'), emit: gtf

  script:
  """
  wget ${params.gtf} | gunzip > genes.gtf
  bash ${baseDir}/bin/prep_resources/filter_gtf.sh genes.gtf filtered_genes.gtf
  """

  stub:
  """
  touch filtered_genes.gtf
  """

}
