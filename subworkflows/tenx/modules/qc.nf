// QC module

nextflow.enable.dsl = 2

//

process QC {

    tag "${sample_name}" 
    label 'MiTo'
    publishDir "${params.output_folder}/${sample_name}", mode: 'copy'   
    
    input:
    tuple val(sample_name), path(filtered)  
    
    output:
    tuple val(sample_name), path('adata.h5ad'), emit: gene_expression_matrix  
    tuple val(sample_name), path('cell_barcodes.txt'), emit: cell_barcodes    
    
    script:

    // Handle CLI from params (only conditional for null-defaulting parameters)
    def path_meta = params.path_meta ? "--path_meta ${params.path_meta}" : ""

    """
    python ${baseDir}/bin/preprocess/QC.py \
    --input ${filtered} \
    --sample ${sample_name} \
    --min_nUMIs ${params.min_nUMIs} \
    --min_n_genes ${params.min_n_genes} \
    --max_perc_mt ${params.max_perc_mt} \
    ${path_meta}
    """

    stub: 
    """
    touch adata.h5ad
    touch cell_barcodes.txt
    """

}
