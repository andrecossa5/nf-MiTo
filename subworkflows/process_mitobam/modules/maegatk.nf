// MAEGATK module

nextflow.enable.dsl = 2

//

process MAEGATK {

  tag "${sample_name}: ${batch_bam.simpleName}"
  label 'maegatk'

  input:
  tuple val(sample_name), path(batch_bam), path(batch_cbs)
  path(reference)

  output:
  tuple val(sample_name),
    val(batch_bam.simpleName),
    path("${batch_bam.simpleName}.A.txt"),
    path("${batch_bam.simpleName}.C.txt"),
    path("${batch_bam.simpleName}.G.txt"),
    path("${batch_bam.simpleName}.T.txt"),
    path("${batch_bam.simpleName}.coverage.txt"),
    path("${batch_bam.simpleName}.depth.txt"), emit: tables

  script:
  def flags = params.maegatk_proper_paired ? '--proper_paired' : ''
  """
  maegatk_batch.sh ${batch_bam} ${batch_bam.simpleName} ${reference[0]} ${params.barcode_tag} \\
    --umi_tag ${params.maegatk_umi_tag} \\
    --NHmax ${params.maegatk_NHmax} \\
    --NMmax ${params.maegatk_NMmax} \\
    --min_reads ${params.maegatk_min_reads} \\
    --base_qual ${params.maegatk_base_qual} \\
    --alignment_quality ${params.maegatk_alignment_quality} \\
    ${flags}
  """

  stub:
  """
  for l in A C G T coverage depth; do touch ${batch_bam.simpleName}.\$l.txt; done
  """

}

//

process COLLAPSE_MAEGATK {

  tag "${sample_name}"
  label 'maegatk'
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

  output:
  tuple val(sample_name), path('tables/*'), emit: ch_output
 
  script:
  """ 
  cat *.A.txt | gzip --fast > A.txt.gz
  cat *.C.txt | gzip --fast > C.txt.gz
  cat *.T.txt | gzip --fast > T.txt.gz
  cat *.G.txt | gzip --fast > G.txt.gz
  cat *.coverage.txt | gzip --fast > coverage.txt.gz
  cat *.depth.txt | gzip --fast > depth.txt.gz
  mkdir tables 
  mv *.gz tables
  """

  stub:
  """
  mkdir tables
  for l in A C G T coverage depth; do touch tables/\$l.txt.gz; done
  """

}
