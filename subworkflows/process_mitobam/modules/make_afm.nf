// MAKE_AFM module

nextflow.enable.dsl = 2

//

process MAKE_AFM {

  tag "${sample_name}" 
  label 'MiTo'
  publishDir "${params.output_folder}/${sample_name}", mode: 'copy'

  input:
  tuple val(sample_name), path(path_ch_matrix), path(fasta)

  output:
  tuple val(sample_name), path('afm_unfiltered.h5ad'), emit: afm 

  script:
  """
  python ${baseDir}/bin/preprocess/make_AF_matrix.py \
  --path_ch_matrix ${path_ch_matrix} \
  --pp_method ${params.pp_method} \
  --scLT_system ${params.scLT_system} \
  --path_meta ${params.path_meta} \
  --sample ${sample_name} 
  """

  stub: 
  """
  touch afm_unfiltered.h5ad
  """

}
