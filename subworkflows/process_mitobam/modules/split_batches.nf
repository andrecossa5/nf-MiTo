// SPLIT_BATCHES module: library MT BAM -> one BAM per batch of cells, in one pass

nextflow.enable.dsl = 2

//

process SPLIT_BATCHES {

  tag "${sample_name}"
  label 'maegatk'

  input:
  tuple val(sample_name), path(mitobam), path(cell_barcodes)

  output:
  tuple val(sample_name), path('batch_*.bam'), path('batch_*.csv'), emit: batches

  script:
  """
  split_bam_by_cb.py ${mitobam} ${cell_barcodes} ${params.batch_size} --barcode_tag ${params.barcode_tag}
  """

  stub:
  """
  for i in 0000 0001 0002; do touch batch_\$i.bam; echo CELL\$i > batch_\$i.csv; done
  """

}
