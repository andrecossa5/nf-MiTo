// PREP_GENOME module

nextflow.enable.dsl = 2

//

process PREP_GENOME {

  tag "${sample_name}"
  label 'scLT'

  output:
  path('masked_genome.fa'), emit: genome

  script:
  """
  wget ${params.reference_genome} | gunzip > genome.fa
  wget ${params.NUMTs_regions} -O hg38.full.blacklist.bed
  bedtools maskfasta -fi genome.fa -bed hg38.full.blacklist.bed -fo masked_genome.fa
  """

  stub:
  """
  touch masked_genome.fa
  """

}
