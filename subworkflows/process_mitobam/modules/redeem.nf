// REDEEMV modules (RedeemV consensus calling, per chunk of cell barcodes)
nextflow.enable.dsl = 2

//

process REDEEM_CONSENSUS {

  tag "${sample_name}: ${batch_bam.simpleName}"
  label 'redeemv'

  input:
  tuple val(sample_name), path(batch_bam), path(batch_cbs)
  tuple path(fasta), path(fai)

  output:
  tuple val(sample_name), path("${batch_bam.simpleName}.RawGenotypes.*.txt"), path("${batch_bam.simpleName}.QualifiedTotalCts"), emit: tables

  script:
  """
  redeem_consensus.py ${batch_bam} ${batch_cbs} ${fasta} ${batch_bam.simpleName} \\
    --barcode_tag ${params.barcode_tag} \\
    --base_quality ${params.redeem_base_quality} \\
    --threads ${task.cpus}
  """

  stub:
  """
  for l in Total VerySensitive Sensitive Specific; do touch ${batch_bam.simpleName}.RawGenotypes.\$l.txt; done
  touch ${batch_bam.simpleName}.QualifiedTotalCts
  """

}

//

process FINALIZE_REDEEM {

  tag "${sample_name}"
  label 'redeemv'
  publishDir "${params.output_folder}", mode: 'copy', saveAs: { fn -> file(fn).name }

  input:
  tuple val(sample_name), path(genotypes), path(counts)

  output:
  tuple val(sample_name), path('tables/*'), emit: ch_output

  script:
  """
  redeem_finalize.py tables \\
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
