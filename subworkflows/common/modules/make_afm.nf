// MAKE_AFM module

nextflow.enable.dsl = 2

//

process MAKE_AFM {

  tag "${sample_name}" 
  publishDir "${params.output_folder}/${sample_name}", mode: 'copy'

  input:
  tuple val(sample_name), path(path_ch_matrix), path(fasta)

  output:
  tuple val(sample_name), path('afm_unfiltered.h5ad'), emit: afm 

  script:
  """
  python ${baseDir}/bin/raw_reads/make_AF_matrix.py \
  # ...
  """

  stub: 
  """
  touch afm_unfiltered.h5ad
  """

}
