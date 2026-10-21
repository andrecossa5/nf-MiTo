// scATAC modules: lanes -> barcoded reads -> bwa-mem2 -> MT BAM

nextflow.enable.dsl = 2

//

process MERGE_ATAC_LANES {

  tag "${sample_name}"
  label 'mgatk'

  input:
  tuple val(sample_name), path(in_folder)

  output:
  tuple val(sample_name), path("R1.fastq.gz"), path("R2.fastq.gz"), path("R3.fastq.gz"), emit: reads

  script:
  """
  # cellranger-atac names: R1 genomic, R2 cell barcode (16 nt), R3 genomic. Lanes are concatenated in the same order.
  for r in R1 R2 R3; do
    ls ${in_folder}/*_\${r}_*.fastq.gz > \${r}.list
    [ -s \${r}.list ] || { echo "No *_\${r}_*.fastq.gz in ${in_folder}" >&2; exit 1; }
    cat \$(cat \${r}.list) > \${r}.fastq.gz
  done
  [ \$(wc -l < R1.list) -eq \$(wc -l < R2.list) ] && [ \$(wc -l < R1.list) -eq \$(wc -l < R3.list) ] \\
    || { echo "Different numbers of R1/R2/R3 files" >&2; exit 1; }
  """

  stub:
  """
  touch R1.fastq.gz R2.fastq.gz R3.fastq.gz
  """

}

//

process ATAC_ADD_BARCODE {

  tag "${sample_name}"
  label 'mgatk'

  input:
  tuple val(sample_name), path(R1), path(R2), path(R3), path(cell_barcodes)

  output:
  tuple val(sample_name), path("barcoded_R1.fastq.gz"), path("barcoded_R3.fastq.gz"), emit: reads

  script:
  """
  atac_add_barcode.py ${R1} ${R2} ${R3} ${cell_barcodes} barcoded --rule ${params.atac_barcode_rule} --compress_threads ${Math.max(1, task.cpus - 2)}
  """

  stub:
  """
  touch barcoded_R1.fastq.gz barcoded_R3.fastq.gz
  """

}

//

process FETCH_MULTIOME_WHITELISTS {

  label 'mgatk'

  output:
  tuple path('WhiteList.RNA'), path('WhiteList.ATAC'), emit: whitelists

  script:
  """
  wget -q -O WhiteList.RNA ${params.multiome_whitelist_rna}
  wget -q -O WhiteList.ATAC ${params.multiome_whitelist_atac}
  """

  stub:
  """
  touch WhiteList.RNA WhiteList.ATAC
  """

}

//

process RNA2ATAC_BARCODES {

  tag "${sample_name}"
  label 'mgatk'

  input:
  tuple val(sample_name), path(rna_cells)
  tuple path(wl_rna), path(wl_atac)

  output:
  tuple val(sample_name), path("atac_cells.tsv"), emit: cell_barcodes

  script:
  """
  rna2atac_barcodes.py ${rna_cells} ${wl_rna} ${wl_atac} atac_cells.tsv
  """

  stub:
  """
  touch atac_cells.tsv
  """

}

//

process BWAMEM2_INDEX {

  label 'mgatk'
  publishDir "${params.output_folder}/resources", mode: 'copy'

  input:
  path(genome)

  output:
  path('bwamem2_index'), emit: index

  script:
  """
  mkdir bwamem2_index
  bwa-mem2 index -p bwamem2_index/genome ${genome}
  """

  stub:
  """
  mkdir bwamem2_index
  """

}

//

process BWAMEM2_MT {

  tag "${sample_name}"
  label 'mgatk'

  input:
  tuple val(sample_name), path(R1), path(R3)
  path(index)

  output:
  tuple val(sample_name), path("mitobam.bam"), emit: bam

  script:
  """
  bwa-mem2 mem -t ${task.cpus} ${index}/genome ${R1} ${R3} \\
    | atac_tag_bam.py ${params.string_MT} ${params.atac_min_mapq} ${params.barcode_tag} unsorted.bam
  samtools sort -@ ${task.cpus} -o mitobam.bam unsorted.bam
  samtools index -@ ${task.cpus} mitobam.bam
  rm unsorted.bam
  """

  stub:
  """
  touch mitobam.bam
  """

}
