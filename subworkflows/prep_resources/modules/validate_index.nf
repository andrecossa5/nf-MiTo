// VALIDATE_INDEX module

nextflow.enable.dsl = 2

//

process VALIDATE_INDEX {

  label 'scLT'

  output:
  path('STAR_index'), emit: star_index
  path('genome.fa'), emit: genome

  script:
  """
  python ${baseDir}/bin/prep_resources/validate_index.py ${params.prebuilt_STAR_index}
  ln -s ${params.prebuilt_STAR_index} STAR_index
  ln -s ${params.prebuilt_STAR_index}/Genome genome.fa
  """

  stub:
  """
  mkdir STAR_index
  touch genome.fa
  """

}
