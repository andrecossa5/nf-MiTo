// REDEEMV modules (RedeemV consensus calling, per chunk of cell barcodes)
nextflow.enable.dsl = 2

//

process REDEEM_CONSENSUS {

  tag "${sample_name}: ${CBs.simpleName}"
  label 'redeemv'

  input:
  tuple val(sample_name), path(bam), path(CBs)
  tuple path(fasta), path(fai)

  output:
  tuple val(sample_name), path("chunk_*.txt"), path("chunk_*.QualifiedTotalCts"), emit: tables

  script:
  def idx = CBs.simpleName.split('_')[1]
  """
  python ${baseDir}/bin/preprocess/redeem_consensus.py ${bam} ${CBs} ${fasta} chunk_${idx} \\
    --redeemv_home ${params.redeemv_home} \\
    --barcode_tag ${params.barcode_tag} \\
    --base_quality ${params.redeem_base_quality}
  """

  stub:
  def idx = CBs.simpleName.split('_')[1]
  """
  for l in Total VerySensitive Sensitive Specific; do touch chunk_${idx}.RawGenotypes.\$l.txt; done
  touch chunk_${idx}.QualifiedTotalCts
  """

}

//

process FINALIZE_REDEEM {

  tag "${sample_name}"
  label 'redeemv'
  publishDir "${params.output_folder}/${sample_name}/MT.preprocess.ouput", mode: 'copy'

  input:
  tuple val(sample_name), path(genotypes), path(counts)

  output:
  tuple val(sample_name), path('tables'), emit: ch_output

  script:
  """
  python ${baseDir}/bin/preprocess/redeem_finalize.py tables \\
    --p_threshold ${params.redeem_strand_bias_p} \\
    --min_log2fc ${params.redeem_strand_bias_log2fc} \\
    --correction ${params.redeem_strand_bias_correction}
  """

  stub:
  """
  mkdir tables
  for l in Total VerySensitive Sensitive Specific; do touch tables/RawGenotypes.\$l.StrandBalance.gz; done
  touch tables/QualifiedTotalCts.gz tables/StrandBiasBlackList
  """

}
