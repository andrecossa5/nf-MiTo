# nf-MiTo Parameter Configuration Guide

This guide provides detailed explanations and best practices for configuring nf-MiTo parameters for different experimental scenarios.

## Quick Reference

| Scenario | Parameter File | Description |
|----------|---------------|-------------|
| Standard MAESTER | `example_maester_basic.json` | Default settings for MAESTER experiments |
| Low Coverage | `example_high_sensitivity.json` | Optimized for low-depth data |
| High Quality | `example_high_stringency.json` | Stringent filtering for high-quality data |
| Parameter Search | `example_parameter_tuning.json` | Grid search configuration |
| Cas9 System | `example_cas9.json` | CRISPR/Cas9 specific settings |
| scWGS | `example_scwgs.json` | Single-cell WGS configuration |
| Benchmarking | `example_benchmarking.json` | Method comparison setup |

## Parameter Categories

### 1. Input/Output Parameters

#### `raw_data_input_type`
- **fastq**: Raw FASTQ files requiring full preprocessing
- **fastq, MAESTER**: MAESTER-specific FASTQ preprocessing
- **mitobam**: Pre-aligned mitochondrial BAM files

#### `output_folder`
- Always required
- Use absolute paths
- Ensure sufficient disk space

### 2. Quality Control Parameters

#### Cell Filtering Thresholds

| Parameter | Low Quality | Standard | High Quality | Description |
|-----------|-------------|----------|--------------|-------------|
| `min_nUMIs` | 200 | 500 | 1000 | Minimum UMIs per cell |
| `min_n_genes` | 100 | 250 | 500 | Minimum genes per cell |
| `max_perc_mt` | 0.20 | 0.15 | 0.10 | Maximum mitochondrial % |
| `n_mads` | 4 | 3 | 2 | MADs for outlier detection |

**Guidelines:**
- **Low quality data**: Relax thresholds to retain more cells
- **High quality data**: Stricter thresholds for cleaner analysis
- **Standard settings**: Good starting point for most datasets

### 3. Variant Detection Parameters

These are the most critical parameters for variant detection quality:

#### Core Detection Parameters

| Parameter | Sensitive | Standard | Stringent | Description |
|-----------|-----------|----------|-----------|-------------|
| `min_n_positive` | 3 | 5 | 10 | Min cells with variant |
| `af_confident_detection` | 0.01 | 0.02 | 0.05 | AF threshold for confidence |
| `min_n_confidently_detected` | 2 | 2 | 5 | Min confident detections |
| `min_mean_AD_in_positives` | 1.0 | 1.25 | 2.0 | Mean allelic depth in positives |
| `t_prob` | 0.6 | 0.7 | 0.8 | Probability threshold |
| `min_AD` | 1 | 2 | 3 | Minimum allelic depth |
| `min_cell_prevalence` | 0.02 | 0.05 | 0.1 | Minimum cell prevalence |

#### Coverage and Quality

| Parameter | Permissive | Standard | Strict | Description |
|-----------|------------|----------|--------|-------------|
| `min_cov` | 3 | 5 | 10 | Minimum coverage |
| `min_var_quality` | 20 | 30 | 40 | Minimum variant quality |
| `min_cell_number` | 3 | 5 | 10 | Min cells with variant |

### 4. System-Specific Recommendations

#### MAESTER
```json
{
    "scLT_system": "MAESTER",
    "pp_method": "maegatk",
    "af_confident_detection": 0.02,
    "min_mean_AD_in_positives": 1.25,
    "distance_metric": "weighted_jaccard"
}
```

#### Cas9
```json
{
    "scLT_system": "Cas9",
    "af_confident_detection": 0.03,
    "min_mean_AD_in_positives": 1.5,
    "distance_metric": "hamming",
    "fgbio_base_error_rate_mito": 0.3
}
```

#### scWGS
```json
{
    "scLT_system": "scWGS",
    "spatial_metrics": true,
    "filter_moran": true,
    "min_cov": 8,
    "tree_algorithm": "iqtree"
}
```

### 5. Phylogeny Parameters

#### Distance Metrics
- **weighted_jaccard**: Best for most MT-SNV data (default)
- **jaccard**: Binary similarity, good for sparse data
- **hamming**: Good for Cas9/editing systems
- **cosine**: Alternative for continuous data

#### Tree Algorithms
- **cassiopeia**: Fast, good for large datasets (default)
- **iqtree**: Maximum likelihood, more accurate
- **mpboot**: Fast bootstrapping

#### Bootstrap Settings
- **Quick analysis**: 100 replicates
- **Standard**: 500 replicates  
- **Publication**: 1000+ replicates

### 6. Parameter Tuning Strategy

#### Recommended Tuning Parameters
Focus tuning on these high-impact parameters:

```json
{
    "min_n_positive": [3, 5, 7, 10],
    "af_confident_detection": [0.01, 0.02, 0.03, 0.05],
    "t_prob": [0.6, 0.7, 0.8],
    "bin_method": ["MiTo", "vanilla"]
}
```

#### Tuning Strategy
1. **Start broad**: Test wide parameter ranges
2. **Narrow down**: Focus on promising regions
3. **Validate**: Test optimal parameters on held-out data

## Usage Examples

### Example 1: Basic MAESTER Analysis
```bash
nextflow run main.nf \
    -profile docker \
    -params-file params/examples/example_maester_basic.json \
    --raw_data_input samples.csv \
    --ref /path/to/reference
```

### Example 2: Parameter Tuning
```bash
nextflow run main.nf \
    -entry TUNE \
    -profile docker \
    -params-file params/examples/example_parameter_tuning.json \
    --afm_input afm_samples.csv
```

### Example 3: High-Sensitivity Analysis
```bash
nextflow run main.nf \
    -entry INFER \
    -profile docker \
    -params-file params/examples/example_high_sensitivity.json \
    --afm_input afm_samples.csv
```

## Troubleshooting Parameter Issues

### Too Few Variants Detected
- Decrease `min_n_positive`
- Decrease `af_confident_detection`
- Decrease `t_prob`
- Decrease `min_cell_prevalence`

### Too Many Noisy Variants
- Increase `min_n_positive`
- Increase `af_confident_detection`
- Increase `min_var_quality`
- Increase `min_mean_AD_in_positives`

### Poor Tree Quality
- Increase `n_boot_replicates`
- Try different `distance_metric`
- Increase variant filtering stringency
- Check for batch effects

### Memory Issues
- Reduce `CBs_chunk_size`
- Use fewer bootstrap replicates
- Consider splitting large datasets

## Parameter Validation

Before running the full pipeline:

1. **Check parameter combinations**: Ensure compatible settings
2. **Validate file paths**: Confirm all input files exist
3. **Test with subset**: Run on small dataset first
4. **Monitor resources**: Check memory and compute requirements

## Best Practices

1. **Document parameters**: Save configuration files with descriptive names
2. **Version control**: Track parameter changes with your data
3. **Systematic testing**: Use TUNE workflow for new datasets
4. **Validate results**: Check output quality metrics
5. **Reproducibility**: Use exact parameter files for publication

## Advanced Configuration

### Custom Parameter Combinations

For complex experiments, you can combine aspects of different configurations:

```json
{
    "_comment": "Custom high-sensitivity Cas9 analysis",
    "scLT_system": "Cas9",
    "min_n_positive": 3,
    "af_confident_detection": 0.015,
    "t_prob": 0.65,
    "distance_metric": "hamming",
    "n_boot_replicates": 500
}
```

### Conditional Parameters

Some parameters only apply to specific workflows:
- `path_tuning`: Only for INFER after TUNE
- `lineage_column`: Only if metadata available
- `spatial_metrics`: Only for spatial datasets

## Getting Help

If you're unsure about parameter settings:

1. **Start with examples**: Use provided example configurations
2. **Run TUNE workflow**: Let the pipeline optimize parameters
3. **Check documentation**: Review parameter descriptions in schema
4. **Community support**: Ask questions via GitHub Issues