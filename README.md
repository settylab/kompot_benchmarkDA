# BenchmarkDA: Differential Abundance Testing Framework

A modular framework for benchmarking differential abundance (DA) methods on single-cell data.

## Quick Start

```bash
# 1. Setup environment
bash setup_environment.sh --minimal

# 2. Run complete pipeline
./cli.sh

# 3. Check status
./cli.sh status
```

## Methods

**Python Methods**: `meld`, `mellon`, `kompot` (multiple variants available)

**R Methods**: `milo`, `daseq`, `cydar`, `louvain`

## Datasets

**Synthetic**: `linear`, `branch`, `cluster` - Generated automatically

**Real**: `covid19-pbmc`, `bcr-xl`, `levine32`, `pancreas` - [Download here](https://drive.google.com/drive/folders/15wWFD5FMe0VdzN1pUnaUUpQ17OXkeebH)

## Common Commands

```bash
# Run specific dataset and method
./cli.sh --datasets linear --methods meld benchmark

# Generate labels with filters
./cli.sh --datasets linear --populations M1,M2 --seeds 43 labels

# Generate labels via SLURM array (one job per population)
./cli.sh --datasets bcr-xl labels --slurm

# Preprocess datasets
./cli.sh --datasets linear,branch preprocess

# Submit benchmarks to SLURM
./cli.sh --datasets linear --methods python benchmark --slurm

# Dry run (see what would execute)
./cli.sh --datasets linear --methods python --dry-run benchmark
```

## Architecture

```
benchmarkDA_private/
├── cli.sh                    # Main interface
├── config/
│   ├── dataset_config.py     # Dataset parameters
│   └── method_config.py      # Method configurations
├── bin/                      # Executables
│   ├── generate_labels.py    # Label generation
│   ├── direct_benchmark.py   # Method execution
│   └── run_method.py         # Unified method wrapper
├── lib/                      # Shared utilities
├── methods/                  # Method implementations
│   ├── meld/
│   ├── mellon/
│   ├── kompot/
│   └── ...
├── data/                     # Input data
│   ├── synthetic/
│   └── real/
└── benchmark/                # Results
```

## Pipeline Steps

1. **Preprocess**: Generate PCA embeddings → compute DM from PCA
2. **Labels**: Generate synthetic condition labels with batch effects
3. **Benchmark**: Run DA methods on all combinations

The CLI handles environment activation automatically.

## Configuration

Edit `config/dataset_config.py` to add datasets:

```python
DATASET_CONFIGS = {
    "linear": {
        "pops": ["M1", "M2", "M3"],        # Populations to test
        "batch_vec": [0, 0.75, 1.5],       # Batch effect levels
        "pop_col": "celltype",              # Population column name
        "n_dm": 10                          # DM components
    }
}
```

Edit `config/method_config.py` to configure methods.

## Results

Results saved to:
```
benchmark/{synthetic|real}/{dataset}/{dataset}-{pop}-{enr}-{seed}-{batch}-{balance}/
```

Each directory contains:
- `*.DAresults.{method}.csv` - Method-specific results
- Metadata and supporting files

## Status Checking

```bash
# Overall status
./cli.sh status

# Specific dataset
./cli.sh --datasets linear status
```

## Filtering Options

Apply filters to labels and benchmarks:

```bash
--populations M1,M2        # Specific populations
--seeds 43,44              # Specific seeds
--enrichments 0.75,0.95    # Enrichment levels
--batch-sds 0,0.75         # Batch standard deviations
```

Example:
```bash
./cli.sh --datasets linear \
         --populations M1,M2 \
         --seeds 43 \
         --enrichments 0.75 \
         --batch-sds 0,0.75 \
         labels benchmark
```

## SLURM Execution

### Label Generation Array Jobs

Submit label generation as SLURM array jobs (one task per population):

```bash
# All populations for a dataset (creates array job)
./cli.sh --datasets bcr-xl labels --slurm

# Multiple datasets (creates one array job per dataset)
./cli.sh --datasets linear,branch labels --slurm

# With filters (only generates labels for specified populations)
./cli.sh --datasets bcr-xl --populations CD4_T-cells,CD8_T-cells labels --slurm

# Custom SLURM options
./cli.sh labels --slurm --sbatch-options "--partition=largenode --mem=32G"
```

Default SLURM settings for labels:
- `--cpus-per-task=8`
- `--time=2-00:00:00`
- Array splits by population automatically

**Note**: All filtering options (`--seeds`, `--enrichments`, `--batch-sds`, `--only-missing`) are automatically forwarded to each array task.

### Benchmark Execution

```bash
# Submit all methods
./cli.sh --datasets linear benchmark --slurm

# Specific methods
./cli.sh --methods milo benchmark --slurm

# Custom SLURM options
./cli.sh benchmark --slurm --sbatch-options "--partition=largenode --mem=64G"
```

## Environment

The `benchmarkda` environment is detected and activated automatically by `cli.sh`.

## Troubleshooting

**Environment not found:**
```bash
bash setup_environment.sh --minimal
```

**Permission denied:**
```bash
chmod +x cli.sh setup_environment.sh
```

**Check available methods:**
```bash
./cli.sh --help
```
