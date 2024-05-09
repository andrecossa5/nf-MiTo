#!/bin/bash -ue
python /Users/IEO5505/Desktop/MI_TO/phylo_inference/bin/filter_variants/get_stats.py     --path_data /Users/IEO5505/Desktop/AML_clonal_reconstruction/data/samples     --sample_name sAML1     --min_site_cov 5     --min_var_quality 30     --min_frac_negative 0.7     --min_n_positive 2     --low_confidence_af 0.001     --high_confidence_af 0.01     --min_prevalence_low_confidence_af 0.01     --min_cells_high_confidence_af 2     --lineage_column malignant_class_occupancy     --solver UPMGA     --metric jaccard     --ncores 4     --path_priors /Users/IEO5505/Desktop/AML_clonal_reconstruction/data/vars_df/priors.csv     --path_meta /Users/IEO5505/Desktop/AML_clonal_reconstruction/data/meta/cells_meta.csv     --job_id 2
