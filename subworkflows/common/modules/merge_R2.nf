// MERGE_R2 module

nextflow.enable.dsl = 2

//

process MERGE_R2 {

  tag "${sample_name}"

  input:
  tuple val(sample_name), path(in_folder), val(library)

  output:
  tuple val(sample_name), path("R2_raw.fastq.gz"), emit: R2

  script:
  """
  zcat ${in_folder}/*R2*.fastq.gz | pigz --fast -p ${task.cpus} > R2_raw.fastq.gz
  """

  stub:
  """
  touch R2_raw.fastq.gz
  """

}
