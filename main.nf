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
    
    DESCRIPTION:
        nf-MiTo is a comprehensive Nextflow pipeline for mitochondrial SNV-based 
        single-cell lineage tracing. It supports multiple lineage tracing systems
        and provides end-to-end analysis from raw data to phylogenetic trees.
    
    USAGE:
        nextflow run main.nf [options]
        nextflow run main.nf -params-file <params.json> [options]
    
    GLOBAL OPTIONS:
        --help                         Show this help message and exit
        --version                      Show version information and exit
        
    ENTRY POINTS (choose one):
        [default]                      Full pipeline: preprocess → infer → annotate
        -entry PREPROCESS              Data preprocessing only (fastq/bam → AFM)
        -entry TUNE                    Parameter optimization (grid search)
        -entry EXPLORE                 Exploratory analysis and visualization
        -entry INFER                   Lineage inference only (AFM → trees)
        -entry BENCH                   Benchmarking and method comparison
        
    EXECUTION PROFILES (choose one or more):
        -profile docker                Use Docker containers (recommended)
        -profile singularity           Use Singularity containers (HPC)
        -profile conda                 Use Conda environments
        -profile local                 Run locally (dependencies required)
    
    ==========================================================================
    INPUT/OUTPUT OPTIONS:
    ==========================================================================

        --raw_data_input_type          Input data format:
                                         • fastq: Raw FASTQ files
                                         • "fastq, MAESTER": MAESTER FASTQ files  
                                         • mitobam: Pre-aligned mitochondrial BAMs
                                       (default: ${params.raw_data_input_type})
                                       
        --raw_data_input               CSV file with raw data paths
                                       Required for: PREPROCESS entry point
                                       Format depends on raw_data_input_type
                                       
        --afm_input                    CSV file with AFM paths (job_id,sample,afm)
                                       Required for: TUNE, EXPLORE, INFER, BENCH
                                       
        --output_folder                Output directory path (REQUIRED)
                                       Will be created if it doesn't exist
                                       
        --path_meta                    Cell metadata file (optional)
                                       (default: ${params.path_meta})
                                       
        --path_tuning                  Tuning results file from TUNE workflow
                                       Use after running TUNE to apply optimal params
                                       (default: ${params.path_tuning})
    
    ==========================================================================
    REFERENCE GENOME OPTIONS:
    ==========================================================================

        --ref                          Reference genome directory
                                       Required for: PREPROCESS entry point
                                       Should contain STAR index files
                                       
        --string_MT                    Mitochondrial chromosome name in reference
                                       Common values: chrM, MT, chrMT
                                       (default: ${params.string_MT})
                                       
        --whitelist                    10x cell barcode whitelist file
                                       Required for: 10x data preprocessing
                                       One barcode per line (text file)
    
    ==========================================================================
    LINEAGE TRACING SYSTEM:
    ==========================================================================

        --scLT_system                  Technology used for lineage tracing:
                                         • MAESTER: mitochondrial RNA editing (default)
                                         • RedeeM: mitochondrial base editing
                                         • Cas9: CRISPR/Cas9 mitochondrial editing
                                         • scWGS: single-cell whole genome sequencing
                                       (default: ${params.scLT_system})
                                       
        --pp_method                    Preprocessing method for variant calling
                                       (default: ${params.pp_method})
    
    ==========================================================================
    KEY PARAMETERS FOR VARIANT DETECTION (use TUNE to optimize):
    ==========================================================================

        --min_n_positive               Minimum cells with variant for inclusion
                                       Lower = more sensitive, Higher = more specific
                                       (default: ${params.min_n_positive})
                                       
        --af_confident_detection       Allele frequency threshold for confident detection
                                       Range: 0.01-0.05, Lower = more sensitive
                                       (default: ${params.af_confident_detection})
                                       
        --t_prob                       Probability threshold for variant calling
                                       Range: 0.6-0.8, Higher = more stringent
                                       (default: ${params.t_prob})
                                       
        --bin_method                   Binarization method [MiTo, vanilla]
                                       MiTo recommended for most analyses
                                       (default: ${params.bin_method})
    
    ==========================================================================
    PHYLOGENY RECONSTRUCTION:
    ==========================================================================

        --distance_metric              Distance metric for tree building:
                                         • weighted_jaccard: best for MT-SNV data
                                         • jaccard: binary similarity
                                         • hamming: good for editing systems  
                                         • cosine: alternative metric
                                       (default: ${params.distance_metric})
                                       
        --tree_algorithm               Tree reconstruction algorithm:
                                         • cassiopeia: fast, scalable (default)
                                         • iqtree: maximum likelihood
                                         • mpboot: fast bootstrapping
                                       (default: ${params.tree_algorithm})
                                       
        --n_boot_replicates            Bootstrap replicates for support values
                                       More replicates = more accurate support
                                       (default: ${params.n_boot_replicates})
                                       
        --lineage_column               Metadata column for lineage annotation
                                       Used for tree annotation if available
                                       (default: ${params.lineage_column})
    
    ==========================================================================
    CELL AND VARIANT FILTERING:
    ==========================================================================

        --min_nUMIs                    Minimum UMIs per cell (default: ${params.min_nUMIs})
        --min_n_genes                  Minimum genes per cell (default: ${params.min_n_genes})
        --max_perc_mt                  Maximum mitochondrial read % (default: ${params.max_perc_mt})
        --min_cov                      Minimum coverage per variant (default: ${params.min_cov})
        --min_var_quality              Minimum variant quality score (default: ${params.min_var_quality})
        --min_cell_number              Minimum cells per variant (default: ${params.min_cell_number})
    
    ==========================================================================
    EXAMPLES:
    ==========================================================================
    
    # Full pipeline with MAESTER data:
    nextflow run main.nf \\
        -profile docker \\
        --raw_data_input samples.csv \\
        --output_folder results \\
        --ref /path/to/reference \\
        --scLT_system MAESTER
    
    # Parameter tuning (recommended first step):
    nextflow run main.nf \\
        -entry TUNE \\
        -profile docker \\
        --afm_input afm_jobs.csv \\
        --output_folder tune_results
    
    # Inference with optimized parameters:
    nextflow run main.nf \\
        -entry INFER \\
        -profile docker \\
        --afm_input afm_jobs.csv \\
        --path_tuning tune_results/optimal_params.csv \\
        --output_folder infer_results
    
    # Using example parameter files:
    nextflow run main.nf \\
        -profile docker \\
        -params-file params/examples/example_maester_basic.json \\
        --afm_input afm_jobs.csv
        
    # Benchmarking analysis:
    nextflow run main.nf \\
        -entry BENCH \\
        -profile docker \\
        --afm_input afm_jobs.csv \\
        --output_folder benchmark_results
    
    ==========================================================================
    PARAMETER CONFIGURATION:
    ==========================================================================
    
    Example parameter files are provided in params/examples/:
        • example_maester_basic.json    - Standard MAESTER analysis
        • example_high_sensitivity.json - Low-coverage data
        • example_high_stringency.json  - High-quality data  
        • example_parameter_tuning.json - Grid search config
        • example_cas9.json             - Cas9 system settings
        • example_scwgs.json            - scWGS configuration
        • example_benchmarking.json     - Benchmarking setup
    
    Use: nextflow run main.nf -params-file <example.json> [other options]
    
    For detailed parameter guidance, see: docs/parameter_guide.md
    
    ==========================================================================
    WORKFLOW RECOMMENDATIONS:
    ==========================================================================
    
    1. NEW DATASET:
       Step 1: Run TUNE to optimize parameters
       Step 2: Run INFER with optimal parameters
       Step 3: Use EXPLORE for visualization
    
    2. ESTABLISHED PROTOCOL:
       Use example parameter files or previous optimal settings
       
    3. METHOD COMPARISON:
       Use BENCH entry point to compare approaches
    
    4. QUICK ANALYSIS:
       Use example parameter files with INFER entry point
    
    ==========================================================================
    SUPPORT:
    ==========================================================================
    
    Documentation: https://github.com/andrecossa5/nf-MiTo
    Issues: https://github.com/andrecossa5/nf-MiTo/issues
    Parameter Guide: docs/parameter_guide.md
    
    For all parameters: nextflow run main.nf --help | less
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
