// BUILD_INDEX module

nextflow.enable.dsl = 2

//

process BUILD_INDEX {

  tag "${sample_name}"
  label 'scLT'

  input:
  path(genome)
  path(gtf)

  output:
  path('STAR_index'), emit: STAR_index

  script:
  """
  STAR \
    --runMode genomeGenerate \
    --genomeDir STAR_index \
    --genomeFastaFiles ${genome} \
    --sjdbGTFfile ${gtf} \
    --runThreadN ${task.cpus}
  """

  stub:
  """
  mkdir STAR_index
  """

}
