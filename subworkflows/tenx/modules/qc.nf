// QC module

nextflow.enable.dsl = 2

//

process QC {

    tag "${sample_name}" 
    publishDir "${params.output_folder}/${sample_name}", mode: 'copy'   
    
    input:
    tuple val(sample_name), path(filtered)  
    
    output:
    tuple val(sample_name), path('adata.h5ad'), emit: gene_expression_matrix  
    tuple val(sample_name), path('cell_barcodes.txt'), emit: cell_barcodes    
    
    // Handle CLI:
    def min_nUMIs = params.min_nUMIs ? "--min_nUMIs ${params.min_nUMIs}" : ""
    def min_n_genes = params.min_n_genes ? "--min_n_genes ${params.min_n_genes}" : ""
    def max_perc_mt = params.max_perc_mt ? "--max_perc_mt ${params.max_perc_mt}" : ""
    def path_meta = params.path_meta ? "--path_meta ${params.path_meta}" : ""
    def n_mads = params.n_mads ? "--n_mads ${params.n_mads}" : ""

    script:
    """
    python ${baseDir}/bin/preprocess/QC.py \
    --input ${filtered} \
    --sample ${sample_name} \
    ${min_nUMIs} \
    ${min_n_genes} \
    ${max_perc_mt} \
    ${path_meta} \
    ${n_mads}
    """

    stub: 
    """
    touch adata.h5ad
    touch cell_barcodes.txt
    """

}
