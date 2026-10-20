// MGATK modules (mgatk bcall: one BAM per cell, counted with mgatk's per-sample function)
nextflow.enable.dsl = 2

//

process MGATK {

  tag "${sample_name}: ${batch_bam.simpleName}"
  label 'mgatk'

  input:
  tuple val(sample_name), path(batch_bam), path(batch_cbs)
  tuple path(fasta), path(fai)

  output:
  tuple val(sample_name), val(batch_bam.simpleName),
    path("${batch_bam.simpleName}.A.txt"),
    path("${batch_bam.simpleName}.C.txt"),
    path("${batch_bam.simpleName}.G.txt"),
    path("${batch_bam.simpleName}.T.txt"),
    path("${batch_bam.simpleName}.coverage.txt"),
    path("${batch_bam.simpleName}.depth.txt"), emit: tables

  script:
  def flags = [
    params.mgatk_keep_duplicates ? '--keep_duplicates' : '',
    params.mgatk_handle_overlap ? '--handle_overlap' : '',
    params.mgatk_proper_paired ? '--proper_paired' : ''
  ].findAll { it }.join(' ')
  """
  mgatk_batch.py ${batch_bam} ${batch_bam.simpleName} ${fasta} \\
    --barcode_tag ${params.barcode_tag} \\
    --umi_tag '${params.mgatk_umi_tag}' \\
    --base_qual ${params.mgatk_base_qual} \\
    --alignment_quality ${params.mgatk_alignment_quality} \\
    --NHmax ${params.mgatk_NHmax} \\
    --NMmax ${params.mgatk_NMmax} \\
    ${flags}
  """

  stub:
  """
  for l in A C G T coverage depth; do touch ${batch_bam.simpleName}.\$l.txt; done
  """

}

//

process COLLAPSE_MGATK {

  tag "${sample_name}"
  label 'mgatk'
  publishDir "${params.output_folder}", mode: 'copy', saveAs: { fn -> file(fn).name }

  input:
  tuple val(sample_name),
    val(batches),
    path(A),
    path(C),
    path(G),
    path(T),
    path(coverage),
    path(depth)
  tuple path(fasta), path(fai)

  output:
  tuple val(sample_name), path('tables/*'), emit: ch_output

  script:
  """
  mgatk_collapse.py ${fasta} tables
  """

  stub:
  """
  mkdir tables
  for l in A C G T coverage depth; do touch tables/\$l.txt.gz; done
  touch tables/refAllele.txt
  """

}
