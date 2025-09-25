// FETCH_WHITELIST module

nextflow.enable.dsl = 2

//

process FETCH_WHITELIST {

  label 'scLT'

  output:
  path('10x_v3_whitelist.txt'), emit: whitelist

  script:
  """
  wget ${params.tenx_whitelist} | gunzip > 10x_v3_whitelist.txt
  """

  stub:
  """
  touch 10x_v3_whitelist.txt
  """

}
