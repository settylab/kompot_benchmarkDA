# BenchmarkDA: Differential Abundance Testing Framework

A modular, configuration-driven framework for benchmarking differential abundance (DA) methods on single-cell data with consistent preprocessing.

## Overview

BenchmarkDA provides a systematic approach to evaluate differential abundance methods across multiple datasets and conditions. The framework handles preprocessing, synthetic label generation, method execution, and comprehensive progress tracking through a single command-line interface.

## Quick Start

```bash
# 1. Setup environment (creates 'benchmarkda' environment)
bash setup_environment.sh --minimal

# 2. Activate the environment (adapt command to your package manager)
conda activate benchmarkda        # For conda users
mamba activate benchmarkda        # For mamba users
micromamba activate benchmarkda   # For micromamba users

# 3. Run commands (environment activated)
python bin/status_report.py      # Check pipeline status
./cli.sh                         # Run complete pipeline
./cli.sh --datasets linear preprocess  # Run specific components
./cli.sh --datasets linear --methods kompot benchmark  # Individual methods
```

**IMPORTANT**: All Python commands require the `benchmarkda` environment. Either:
- Activate the environment first (shown above), then use `python` normally
- Use `./cli.sh` which handles environment automatically
- Use direct environment execution: `conda run -n benchmarkda python ...` (or `mamba run`, `micromamba run`)

## Methods and Datasets

### Python Methods
- **kompot**: Kernel-based differential abundance testing
- **mellon**: Manifold-based methods (multiple variants)
- **meld**: Manifold enhancement of latent dimensions
- **cna**: Conditional neighborhood analysis

### R Methods
- **milo**: Neighborhood-based DA testing
- **daseq**: DA region detection
- **cydar**: Hypersphere-based DA testing
- **louvain**: Clustering-based DA testing

### Datasets
- **Synthetic**: `linear`, `branch`, `cluster` - Generated topologies for benchmarking
- **Real**: `covid19-pbmc`, `bcr-xl`, `levine32`, `pancreas` - Real biological datasets

## Architecture

### Core Features
- **Configuration-driven**: All datasets and methods defined in `config/`
- **Single CLI**: The `cli.sh` interface for all operations
- **Status tracking**: Comprehensive progress monitoring
- **Direct execution**: No intermediate script generation
- **Environment detection**: Works with mamba/conda/micromamba
- **Modular pipeline**: Independent preprocessing, labels, and benchmarking stages

### Directory Structure

```
benchmarkDA_private/
├── cli.sh                           # Main CLI interface
├── setup_environment.sh             # Environment setup
│
├── bin/
│   ├── status_report.py             # Progress tracking and reporting
│   ├── direct_benchmark.py          # Direct method execution
│   ├── environment_utils.sh         # Environment detection utilities
│   ├── dataset_preprocessing.sh     # Dataset preprocessing pipeline
│   └── modified_benchmarkda_dm_all.sh # Label generation pipeline
│
├── config/
│   ├── dataset_config.py            # Dataset parameters and populations
│   └── method_config.py             # Method configurations
│
├── python_method/
│   ├── kompot_bm.py, Mellon_bm.py   # Python method implementations
│   ├── meld_bm.py, CNA_bm.py        # Method implementations
│   ├── generate_bm_data.py          # Shared data generation
│   └── *.py                         # Supporting utilities
│
├── scripts/
│   └── run_DA.r                     # R method interface
│
├── data/
│   ├── synthetic/{linear,branch,cluster}/  # Generated datasets
│   └── real/{covid19-pbmc,pancreas,...}/   # Real datasets
│
├── benchmark/
│   ├── synthetic/{dataset}/{job-id}/  # Synthetic results
│   └── real/{dataset}/{job-id}/       # Real dataset results
│
└── environment_minimal.yml          # Environment specification
```

## Usage

### Pipeline Management
```bash
# Activate environment first (choose your package manager)
conda activate benchmarkda  # or mamba/micromamba activate benchmarkda

# Then run commands normally
python bin/status_report.py --detailed --commands  # Check what needs to be run
./cli.sh                                          # Run complete pipeline
./cli.sh --datasets linear,branch preprocess      # Run specific steps
python bin/status_report.py --dataset levine32 --detailed  # Dataset status
```

### Individual Method Testing
```bash
# Test specific methods on specific datasets
./cli.sh --datasets linear --methods kompot benchmark
./cli.sh --datasets branch --methods mellon benchmark

# Run multiple methods
./cli.sh --datasets linear --methods kompot,mellon benchmark

# Use SLURM for larger runs
./cli.sh --datasets linear --methods python benchmark --slurm
```

### Advanced Options
```bash
# Dry run to see what would be executed
./cli.sh --datasets linear --methods python --dry-run

# Skip missing datasets
./cli.sh --skip-missing preprocess

# Custom embedding selection
./cli.sh --embeddings pca benchmark
```

## Environment Setup

### Automatic Setup
```bash
# Interactive setup with environment detection
bash setup_environment.sh

# Minimal setup (non-interactive)
bash setup_environment.sh --minimal
```

### Manual Setup
```bash
# Using your preferred package manager
mamba env create -f environment_minimal.yml
mamba activate benchmarkda
```

### Environment Detection
The system automatically detects and uses:
- MAMBA_EXE or CONDA_EXE environment variables (preferred)
- Fallback to `micromamba` if available
- Works around IT placeholder scripts for institutional setups

## Status Monitoring

```bash
# Activate environment first
conda activate benchmarkda  # or mamba/micromamba activate benchmarkda

# Comprehensive status report
python bin/status_report.py --detailed

# Summary with suggested commands
python bin/status_report.py --commands

# Save status to JSON
python bin/status_report.py --save

# Individual dataset status
python bin/status_report.py --dataset pancreas --detailed
```

**Note**: The `cli.sh` script automatically handles the environment activation, so you can always use it directly without manually activating.

The status report provides:
- **Data file availability**: Which datasets are present
- **Preprocessing status**: DM/PCA embedding completion
- **Label generation**: Synthetic condition combinations completed
- **Benchmark progress**: Individual method completion by embedding type
- **Suggested commands**: Exact CLI commands to run missing components

## Configuration

### Dataset Configuration
Dataset parameters are defined in `config/dataset_config.py`:

```python
DATASET_CONFIGS = {
    "linear": {
        "pops": ["M1", "M2", "M3", "M4", "M5", "M6", "M7"],  # Cell populations to test
        "batch_vec": [0, 0.75, 1, 1.25, 1.5],               # Batch effect levels
        "pop_col": "celltype",                                # Column name for cell types
        "n_dm": 10                                            # Diffusion map components
    }
}

SEEDS = [43, 44, 45]                    # Random seeds for reproducibility
ENRICHMENT_VALUES = [0.75, 0.85, 0.95] # Population enrichment levels
```

### Method Configuration
Method-specific parameters and execution details are defined in `config/method_config.py`.

## Pipeline Workflow

The benchmarking pipeline consists of three independent stages:

1. **Preprocessing**: Generate PCA and DM embeddings for each dataset
2. **Label Generation**: Create synthetic condition labels (independent of embeddings)
3. **Benchmarking**: Run all DA methods on all population/enrichment/seed combinations

### Label Combinations
For each dataset, labels are generated for all combinations of:
- Populations × Seeds × Enrichment values × Batch settings
- Example: 8 populations × 3 seeds × 3 enrichments × 1 batch = 72 combinations

### Benchmark Jobs
Each method runs on all label combinations for both DM and PCA embeddings.

## Results Structure

Results are organized in a hierarchical structure:
```
benchmark/{synthetic|real}/{dataset}/{job-id}/iteration_0/
```

Where `job-id` follows the format: `{dataset}-{population}-{enrichment}-{seed}-{batch}-{balance}-{embedding}`

## Implementation Details

### Label Generation
Labels are created once per dataset, ensuring consistency across method comparisons.

### Direct Method Execution
Methods are executed directly without intermediate script generation, reducing complexity and improving maintainability.

### Environment Compatibility
The framework includes robust environment detection that works with various conda/mamba configurations and handles institutional IT constraints.

### Status Tracking
Comprehensive progress monitoring tracks completion at the dataset, method, and embedding level, enabling resumable execution and clear progress reporting.

## Performance Characteristics

- **Parallel execution**: Methods run independently
- **Resource efficient**: Only generates necessary embeddings
- **Resumable**: Status tracking allows continuing interrupted runs
- **Modular**: Individual components can be run independently

## Data Requirements

**Synthetic datasets** (linear, branch, cluster): Generated automatically

**Real datasets** (covid19-pbmc, bcr-xl, levine32, pancreas): Must be downloaded separately
- Download from: [Google Drive link](https://drive.google.com/drive/folders/15wWFD5FMe0VdzN1pUnaUUpQ17OXkeebH)
- Place in `data/real/{dataset}/` directories

## Environment Usage Summary

**Three ways to run Python commands:**

1. **Environment activation** (recommended):
   ```bash
   conda activate benchmarkda     # (or mamba/micromamba activate benchmarkda)
   python bin/status_report.py   # Use python normally after activation
   ```

2. **CLI script** (automatic environment handling):
   ```bash
   ./cli.sh preprocess           # Handles environment automatically
   ./cli.sh --datasets linear benchmark
   ./cli.sh status
   ```

3. **Direct environment execution** (no activation needed):
   ```bash
   conda run -n benchmarkda python bin/status_report.py    # For conda users
   mamba run -n benchmarkda python bin/status_report.py    # For mamba users
   micromamba run -n benchmarkda python bin/status_report.py  # For micromamba users
   ```

**Never use system python directly** - commands like `python bin/status_report.py` will fail unless you've first activated the `benchmarkda` environment.

## Troubleshooting

### Environment Issues
```bash
# Environment not found
bash setup_environment.sh --minimal

# Package manager not detected - set environment variables:
export MAMBA_EXE=/path/to/mamba
# or
export CONDA_EXE=/path/to/conda
```

### Permission Issues
```bash
# If you see "permission denied" errors:
chmod +x cli.sh
chmod +x setup_environment.sh
```

### Status Report Issues
```bash
# If preprocessing not detected, check file naming:
ls data/synthetic/linear/linear_DM_*.h5ad
ls data/synthetic/linear/linear_PCA_*.h5ad
```

For questions or issues, consult the status report: `python bin/status_report.py --commands` (after environment activation)