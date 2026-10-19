// MGATK modules (mgatk bcall: one BAM per cell, counted with mgatk's per-sample function)
nextflow.enable.dsl = 2

//

process MGATK {

  tag "${sample_name}: ${cell}"
  label 'mgatk'

  input:
  tuple val(sample_name), val(cell), path(bam)
  tuple path(fasta), path(fai)

  output:
  tuple val(sample_name), val(cell),
    path("${cell}.A.txt"),
    path("${cell}.C.txt"),
    path("${cell}.G.txt"),
    path("${cell}.T.txt"),
    path("${cell}.coverage.txt"),
    path("${cell}.depth.txt"), emit: tables

  script:
  def flags = [
    params.mgatk_keep_duplicates ? '--keep_duplicates' : '',
    params.mgatk_handle_overlap ? '--handle_overlap' : '',
    params.mgatk_proper_paired ? '--proper_paired' : ''
  ].findAll { it }.join(' ')
  """
  mgatk_cell.py ${bam} ${cell} ${fasta} \\
    --umi_tag '${params.mgatk_umi_tag}' \\
    --base_qual ${params.mgatk_base_qual} \\
    --alignment_quality ${params.mgatk_alignment_quality} \\
    --NHmax ${params.mgatk_NHmax} \\
    --NMmax ${params.mgatk_NMmax} \\
    ${flags}
  """

  stub:
  """
  for l in A C G T coverage depth; do touch ${cell}.\$l.txt; done
  """

}

//

process COLLAPSE_MGATK {

  tag "${sample_name}"
  label 'mgatk'
  publishDir "${params.output_folder}/${sample_name}/MT.preprocess.ouput", mode: 'copy'

  input:
  tuple val(sample_name),
    val(cells),
    path(A),
    path(C),
    path(G),
    path(T),
    path(coverage),
    path(depth)
  tuple path(fasta), path(fai)

  output:
  tuple val(sample_name), path('tables'), emit: ch_output

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
