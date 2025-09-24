# BenchmarkDA: Differential Abundance Testing Benchmarks

This repository contains a unified framework for benchmarking differential abundance (DA) testing methods in single-cell data.

## Overview

Differential abundance testing aims to identify cell types or states that change in proportion between conditions in single-cell data. This benchmark evaluates multiple methods:

**Python Methods:**
- MELD
- CNA (Conditional Neighbor Analysis)
- Mellon

**R Methods:**
- Milo
- DAseq
- CyDAR

The benchmark evaluates these methods using both real datasets and synthetic topologies, with various parameters including:
- PCA and diffusion map (DM) embeddings
- Different cell populations
- Various enrichment levels
- Different batch effect strengths
- Multiple random seeds

## Requirements

Our implementation requires a Slurm job scheduler since we need to run thousands of parallel jobs.

### Quick Setup

```bash
# Create the environment (choose one):
bash setup_environment.sh --minimal    # Flexible versions
bash setup_environment.sh --complete   # Pinned versions

# Run the benchmark:
bash main.sh
```

### Dependencies

- **Python environment** (for Python methods):
  - Environment name: `benchmarkda` (user-agnostic)
  - Managed with conda/mamba/micromamba
  - Two environment options:
    - `environment_complete.yml`: Pinned versions (recommended for reproducibility)
    - `environment_minimal.yml`: Flexible versions (better cross-platform compatibility)
  - Includes kompot from latest master: `git+https://github.com/settylab/kompot.git@master`

- **R environment** (for R methods):
  - Managed with renv
  - Configuration in `renv.lock`
  - Required packages: miloR, DAseq, cydar, SingleCellExperiment

## Datasets

The benchmark uses two types of datasets:

1. **Synthetic topologies**:
   - Linear
   - Branch
   - Cluster

2. **Real datasets** (must be downloaded separately, e.g., from [Google Drive](https://drive.google.com/drive/folders/15wWFD5FMe0VdzN1pUnaUUpQ17OXkeebH)):
   - bcr-xl
   - covid19-pbmc
   - levine32
   - pancreas

## Running the Benchmark

The entire workflow is orchestrated by the `main.sh` script, which:

1. Sets up required Python and R environments
2. Creates necessary directory structure
3. Preprocesses datasets with PCA and diffusion map embeddings
4. Generates synthetic condition labels with known ground truth
5. Runs all methods across all parameter combinations using Slurm

To run the full benchmark:

```bash
# First-time setup:
bash setup_environment.sh

# Run benchmark:
bash main.sh
```

## Architecture

The benchmark has been refactored for improved reliability and readability:

- **Configuration Files**:
  - `config/dataset_config.py`: Dataset parameters and path templates
  - `config/method_config.py`: Method configurations and command generation

- **Core Components**:
  - `python_method/data_loader.py`: Unified data loading for all methods
  - Python method implementations: `Mellon_bm.py`, `meld_bm.py`, `CNA_bm.py`
  - R method execution: Handled via `scripts/run_DA.r`

- **Execution Framework**:
  - `bin/run_benchmark.py`: Unified script generator for all methods
  - Runtime-generated scripts: Stored in `benchmark_scripts/` directory
  - Slurm job logs: Stored in `SlurmLog/` directory

## Output Structure

Results are organized in the `benchmark/` directory with the following structure:

```
benchmark/
├── dm/                   # Diffusion map results
│   ├── synthetic/
│   │   ├── linear/
│   │   ├── branch/
│   │   └── cluster/
│   └── real/
│       ├── bcr-xl/
│       ├── covid19-pbmc/
│       └── ...
└── pca/                  # PCA results
    ├── synthetic/
    │   ├── linear/
    │   ├── branch/
    │   └── cluster/
    └── real/
        ├── bcr-xl/
        ├── covid19-pbmc/
        └── ...
```

Within each dataset directory, results are further organized by specific parameters:
`[dataset]-[population]-[enrichment]-[seed]-[batch_sd]-[balance]-[embedding]/`

## Developers Guide

To extend the benchmark with new methods:

1. Add method configuration to `config/method_config.py`
2. For Python methods:
   - Implement a wrapper in `python_method/`
   - Use the unified data loading interface
3. For R methods:
   - Add handling to `scripts/run_DA.r`
4. Run the benchmark with your new method
