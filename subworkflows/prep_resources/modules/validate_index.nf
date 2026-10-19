// VALIDATE_INDEX module

nextflow.enable.dsl = 2

//

process VALIDATE_INDEX {

  label 'maegatk'

  output:
  path('STAR_index'), emit: star_index

  script:
  """
  validate_index.py ${params.prebuilt_STAR_index}
  ln -s ${params.prebuilt_STAR_index} STAR_index
  """

  stub:
  """
  mkdir STAR_index
  """

}
