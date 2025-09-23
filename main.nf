// nf-MiTo
nextflow.enable.dsl = 2

// Include subworkflows
include { preprocess } from "./subworkflows/preprocess/main"
include { createPreprocessingChannel } from "./subworkflows/preprocess/main"
include { afm_preprocess } from "./subworkflows/afm_preprocess/main"
include { createAFMChannel } from "./subworkflows/afm_preprocess/main"
include { build_tree } from "./subworkflows/tree_build/main"
include { process_tree } from "./subworkflows/tree_process/main"
include { tune } from "./subworkflows/tune/main"
include { explore } from "./subworkflows/explore/main"
include { benchmark } from "./subworkflows/bench/main"


//


// Help and version
if (params.help) {

    println """\
    ==========================================================================
    nf-MiTo: Integrated pipeline for MT-SNVs-based single-cell lineage tracing
    ==========================================================================
    
    Usage:
        nextflow run main.nf [options]
    
    GLOBAL OPTIONS:
        --help                         Show this help message
        --version                      Show version information
        
    ENTRY POINTS:
        -entry PREPROCESS              Run preprocessing only
        -entry TUNE                    Run parameter tuning
        -entry EXPLORE                 Run exploratory analysis  
        -entry INFER                   Run inference workflow (default)
        -entry BENCH                   Run benchmarking
        
    EXECUTION PROFILES:
        -profile docker                Use Docker containers
        -profile singularity           Use Singularity containers
        -profile conda                 Use Conda environments
        -profile local                 Run locally
    
    ==========================================================================
    INPUT/OUTPUT OPTIONS:
    ==========================================================================

        --raw_data_input_type          Type of raw input data [fastq, "fastq, MAESTER", mitobam] (default: ${params.raw_data_input_type})
        --raw_data_input               Path to CSV file containing raw data input information
        --afm_input                    Path to CSV file containing AFM input information
        --output_folder                Output directory (REQUIRED)
        --path_meta                    Path to metadata file (default: ${params.path_meta})
        --path_tuning                  Path to tuning results file (default: ${params.path_tuning})
    
    ==========================================================================
    REFERENCE GENOME OPTIONS:
    ==========================================================================

        --ref                          Path to reference genome directory
        --string_MT                    Mitochondrial chromosome identifier (default: ${params.string_MT})
        --whitelist                    Path to cell barcode whitelist file
    
    ==========================================================================
    scLT SYSTEM OPTIONS:
    ==========================================================================

        --scLT_system                  Single-cell lineage tracing system [MAESTER, RedeeM, Cas9, scWGS] (default: ${params.scLT_system})
        --pp_method                    Preprocessing method for variant calling (default: ${params.pp_method})
    
    ==========================================================================
    SEQUENCING DATA PREPROCESSING OPTIONS:
    ==========================================================================

        --CBs_chunk_size               Cell barcode chunk size for processing (default: ${params.CBs_chunk_size})
        --fgbio_UMI_consensus_mode     fgbio UMI consensus calling strategy (default: ${params.fgbio_UMI_consensus_mode})
        --fgbio_UMI_consensus_edits    Maximum edit distance for UMI consensus (default: ${params.fgbio_UMI_consensus_edits})
        --fgbio_min_reads_mito         Minimum reads required for mitochondrial consensus (default: ${params.fgbio_min_reads_mito})
        --fgbio_base_error_rate_mito   Base error rate for mitochondrial consensus calling (default: ${params.fgbio_base_error_rate_mito})
        --fgbio_base_quality           Minimum base quality for consensus calling (default: ${params.fgbio_base_quality})
        --fgbio_min_alignment_quality  Minimum alignment quality score (default: ${params.fgbio_min_alignment_quality})
    
    ==========================================================================
    CELL FILTERING OPTIONS:
    ==========================================================================

        --min_nUMIs                    Minimum number of UMIs per cell (default: ${params.min_nUMIs})
        --min_n_genes                  Minimum number of genes per cell (default: ${params.min_n_genes})
        --max_perc_mt                  Maximum percentage of mitochondrial reads (default: ${params.max_perc_mt})
        --n_mads                       Number of median absolute deviations for outlier detection (default: ${params.n_mads})
    
    ==========================================================================
    ALLELE FREQUENCY MATRIX PREPROCESSING OPTIONS:
    ==========================================================================

        --filter_dbs                   Filter doublets/blacklisted variants (default: ${params.filter_dbs})
        --cell_filter                  Cell filtering strategy (default: ${params.cell_filter})
        --filtering                    Variant filtering method (default: ${params.filtering})
        --spatial_metrics              Calculate spatial metrics (default: ${params.spatial_metrics})
        --filter_moran                 Apply Moran's I filtering (default: ${params.filter_moran})
        --min_cell_number              Minimum number of cells with variant (default: ${params.min_cell_number})
        --min_cov                      Minimum coverage (default: ${params.min_cov})
        --min_var_quality              Minimum variant quality score (default: ${params.min_var_quality})
        --min_frac_negative            Minimum fraction of negative cells (default: ${params.min_frac_negative})
        --min_n_positive               Minimum number of positive cells (default: ${params.min_n_positive})
        --af_confident_detection       Allele frequency threshold for confident detection (default: ${params.af_confident_detection})
        --min_n_confidently_detected   Minimum number of confidently detected variants (default: ${params.min_n_confidently_detected})
        --min_mean_AD_in_positives     Minimum mean allelic depth in positive cells (default: ${params.min_mean_AD_in_positives})
        --min_mean_DP_in_positives     Minimum mean depth in positive cells (default: ${params.min_mean_DP_in_positives})
        --t_prob                       Probability threshold (default: ${params.t_prob})
        --min_AD                       Minimum allelic depth (default: ${params.min_AD})
        --min_cell_prevalence          Minimum cell prevalence for variant (default: ${params.min_cell_prevalence})
        --t_vanilla                    Vanilla threshold (default: ${params.t_vanilla})
        --bin_method                   Binarization method (default: ${params.bin_method})
        --k                            Number of neighbors for kNN (default: ${params.k})
        --gamma                        Gamma parameter (default: ${params.gamma})
        --min_n_var                    Minimum number of variants per cell (default: ${params.min_n_var})
    
    ==========================================================================
    PHYLOGENY RECONSTRUCTION OPTIONS:
    ==========================================================================

        --lineage_column               Column name for lineage annotation (default: ${params.lineage_column})
        --K                            Number of clusters (default: ${params.K})
        --distance_metric              Distance metric for phylogenetic analysis [weighted_jaccard, jaccard, hamming, cosine] (default: ${params.distance_metric})
        --tree_algorithm               Tree reconstruction algorithm [cassiopeia, iqtree, mpboot] (default: ${params.tree_algorithm})
        --cassiopeia_solver            Cassiopeia tree solver method [UPMGA, NJ, UPGMA] (default: ${params.cassiopeia_solver})
        --n_boot_replicates            Number of bootstrap replicates (default: ${params.n_boot_replicates})
        --boot_strategy                Bootstrap strategy [feature_resampling, cell_resampling] (default: ${params.boot_strategy})
        --frac_char_resampling         Fraction of characters to resample (default: ${params.frac_char_resampling})
        --support_method               Support calculation method [tbe, felsenstein] (default: ${params.support_method})
        --annotate_tree                Tree annotation method (default: ${params.annotate_tree})
        --max_fraction_unassigned      Maximum fraction of unassigned cells (default: ${params.max_fraction_unassigned})
    
    ==========================================================================
    BENCHMARKING OPTIONS:
    ==========================================================================

        --maxK                         Maximum number of clusters for vireo (default: ${params.maxK})
    
    ==========================================================================
    EXAMPLES:
    ==========================================================================
    
    # Run full pipeline with mitobam input:
    nextflow run main.nf 
        -profile docker,local \\
        --raw_data_input samples.csv \\
        --output_folder results \\
        --ref /path/to/reference 
    
    # Run parameter tuning:
    nextflow run main.nf 
        -entry TUNE 
        -profile docker,local \\
        --afm_input afm_jobs.csv \\
        --output_folder tune_results
    
    # Run inference only:
    nextflow run main.nf \\
        -entry INFER \\
        -profile docker,local \\
        --afm_input afm_jobs.csv \\
        --output_folder infer_results
    
    ==========================================================================
    For more information: https://github.com/andrecossa5/nf-MiTo
    ==========================================================================
    """.stripIndent()
    exit(0)

}

// Version
if (params.version) {
    println "nf-MiTo version: ${workflow.manifest.version}"
    exit(0)
}

// Parameter validation
if (!params.output_folder && !params.help && !params.version) {
    error "Error: --output_folder is required. Use --help for usage information."
}


//


//----------------------------------------------------------------------------//
// nf-MiTo complete entry-points
//----------------------------------------------------------------------------//

// Raw data pre-processing
workflow PREPROCESS {

    ch_preprocessing = createPreprocessingChannel()
    preprocess(ch_preprocessing)

}

//----------------------------------------------------------------------------//

// Lineage inference
workflow TUNE {

    ch_jobs = createAFMChannel()
    tune(ch_jobs)

}

//

workflow EXPLORE {

    ch_jobs = createAFMChannel()
    explore(ch_jobs)

}

//

workflow INFER {

    ch_jobs = createAFMChannel()
    afm_preprocess(ch_jobs)
    build_tree(afm_preprocess.out.input)
    process_tree(afm_preprocess.out.input, build_tree.out.final_tree)

}

//

workflow BENCH {

    ch_jobs = createAFMChannel()
    benchmark(ch_jobs)
    benchmark.out.pickles.view()

}

//


//----------------------------------------------------------------------------//
// nf-MiTo main workflow
//----------------------------------------------------------------------------//

workflow {
    
    println "\n"
    println "This is nf-MiTo, the Nextflow pipeline for MT-SNVs-based scLT."
    println "Usage: nextflow run main.nf -c <config> -params-file <params> -profile <profile> -entry <entry>"
    println "nextflow run main.nf --help for all options and usage."
    println "\n"

    ch_pp = createPreprocessingChannel()
    preprocess(ch_pp)
    afm_preprocess(preprocess.out.afm)
    build_tree(afm_preprocess.out.input)
    process_tree(afm_preprocess.out.input, build_tree.out.final_tree)

}


//
