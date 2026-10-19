// MGATK modules (mgatk tenx counting, run per chunk of cell barcodes)
nextflow.enable.dsl = 2

//

process MGATK {

  tag "${sample_name}: ${CBs.simpleName}"
  label 'mgatk'

  input:
  tuple val(sample_name), path(bam), path(CBs)
  tuple path(fasta), path(fai)

  output:
  tuple val(sample_name), path("chunk_*.txt"), emit: tables

  script:
  def idx = CBs.simpleName.split('_')[1]
  def flags = [
    params.mgatk_keep_duplicates ? '--keep_duplicates' : '',
    params.mgatk_handle_overlap ? '--handle_overlap' : '',
    params.mgatk_proper_paired ? '--proper_paired' : ''
  ].findAll { it }.join(' ')
  """
  python ${baseDir}/bin/preprocess/mgatk_chunk.py ${bam} ${CBs} ${fasta} chunk_${idx} \\
    --barcode_tag ${params.barcode_tag} \\
    --umi_tag '${params.mgatk_umi_tag}' \\
    --base_qual ${params.mgatk_base_qual} \\
    --alignment_quality ${params.mgatk_alignment_quality} \\
    --NHmax ${params.mgatk_NHmax} \\
    --NMmax ${params.mgatk_NMmax} \\
    ${flags}
  """

  stub:
  def idx = CBs.simpleName.split('_')[1]
  """
  for l in A C G T coverage depth; do touch chunk_${idx}.\$l.txt; done
  """

}

//

process COLLAPSE_MGATK {

  tag "${sample_name}"
  label 'mgatk'
  publishDir "${params.output_folder}/${sample_name}/MT.preprocess.ouput", mode: 'copy'

  input:
  tuple val(sample_name), path(chunk_tables)
  tuple path(fasta), path(fai)

  output:
  tuple val(sample_name), path('tables'), emit: ch_output

  script:
  """
  python ${baseDir}/bin/preprocess/mgatk_collapse.py ${fasta} tables
  """

  stub:
  """
  mkdir tables
  for l in A C G T coverage depth; do touch tables/\$l.txt.gz; done
  touch tables/refAllele.txt
  """

}
