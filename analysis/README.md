# BenchmarkDA Results Analysis

Automated analysis framework for differential abundance benchmarking results.

## Features

- **Automatic result discovery** from unified benchmark structure
- **Comprehensive metrics**: AUROC, AUPRC, Pearson/Spearman correlation, MAE, etc.
- **Method comparison** across datasets, populations, and parameters
- **Visualization tools** for performance comparison and ranking
- **Ground truth integration** for accurate metric computation

## Quick Start

```python
from analysis.result_analyzer import ResultAnalyzer, create_default_config

# Create configuration
config = create_default_config(
    datasets=['linear'],
    methods=['meld', 'mellon', 'milo', 'louvain'],
    metrics=['auroc', 'auprc', 'pearson']
)

# Analyze results
analyzer = ResultAnalyzer(config)
metrics_df = analyzer.analyze_all_results()

# View results
print(f"Analyzed {len(metrics_df)} evaluations")
print(metrics_df.groupby('method')['auroc'].mean())
```

## Configuration Options

### Datasets
```python
datasets=['linear', 'branch', 'cluster', 'covid19-pbmc', 'bcr-xl', 'levine32', 'pancreas']
```

### Methods
**Python methods:**
- `meld`, `meld_default` - MELD variants
- `mellon`, `mellon_pca`, `mellon_noSync`, `mellon_corr` - Mellon variants
- `kompot`, `kompot_pca` - Kompot variants

**R methods:**
- `milo` - Milo
- `louvain` - Louvain clustering
- `daseq` - DAseq (requires package installation)
- `cydar` - CyDAR (currently not working)

### Filtering Options
```python
config = create_default_config(
    datasets=['linear'],
    populations=['M1', 'M2'],
    enrichments=[0.75, 0.85, 0.95],
    seeds=[43, 44, 45],
    batch_effects=[0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4],
    methods=['meld', 'mellon', 'milo'],
    metrics=['auroc', 'auprc', 'pearson']
)
```

## Available Metrics

- **Classification**: `auroc`, `auprc`, `auroc_da`, `auroc_neg`, `auroc_pos`
- **Correlation**: `pearson`, `spearman`, `kendall`
- **Error**: `mae`, `mse`, `r2`

## Visualization

### Performance Comparison
```python
fig = analyzer.plot_performance_comparison(
    metrics_df,
    metric='auroc',
    group_by='enrichment',
    figsize=(10, 6)
)
```

### Method Ranking
```python
fig = analyzer.plot_method_ranking(
    metrics_df,
    metric='auroc',
    figsize=(12, 6)
)
```

## Jupyter Notebook

See `notebooks/results_analysis.ipynb` for comprehensive working examples including:
- Result discovery and inspection
- Metrics computation
- Performance visualization
- Method comparison
- Parameter filtering

## Project Structure

```
analysis/
├── result_analyzer.py    # Main analysis framework
└── README.md             # This file

notebooks/
├── results_analysis.ipynb  # Working examples
└── archive/                # Old/deprecated notebooks
```

## Notes

- The analyzer automatically finds the project root by looking for `cli.sh`, `benchmark/`, or `config/dataset_config.py`
- Results must be in the unified structure (no `iteration_0` subdirectory, no embedding suffix)
- Ground truth labels must exist in `data/{synthetic|real}/{dataset}/{job-id}/` for metrics computation
- Method results use standardized `col_0` column for predictions

## Troubleshooting

**No results found:**
- Check that benchmark jobs have completed successfully
- Verify results are in `benchmark/synthetic/` or `benchmark/real/`
- Ensure job directory names match pattern: `{dataset}-{pop}-{enr}-{seed}-{batch}-No`

**Metrics are NaN:**
- Verify ground truth labels exist in `data/` directory
- Check that batch effect values match between results and ground truth
- Ensure method result files have `col_0` column

**Old batch effect results:**
- Previous runs used batch effects: [0, 0.75, 1.0, 1.25, 1.5]
- Current configuration uses: [0, 0.05, 0.1, 0.2, 0.4, 1, 2, 4]
- Filter by batch effects to avoid mixing configurations
