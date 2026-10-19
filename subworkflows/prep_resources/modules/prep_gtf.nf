// PREP_GTF module

nextflow.enable.dsl = 2

//

process PREP_GTF {

  label 'scLT'

  output:
  path('filtered_genes.gtf'), emit: gtf

  script:
  """
  wget -qO- ${params.gtf} | gunzip > genes.gtf
  filter_gtf.sh genes.gtf filtered_genes.gtf
  """

  stub:
  """
  touch filtered_genes.gtf
  """

}
