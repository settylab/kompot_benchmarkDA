# BenchmarkDA Setup Guide

## Quick Start (For Any User)

### 1. Environment Setup

```bash
# Clone the repository (if not already done)
git clone <repository-url>
cd benchmarkDA_private

# Create the environment automatically
bash setup_environment.sh

# Follow the prompts to choose:
#   1) Minimal (flexible versions, better compatibility)
#   2) Complete (pinned versions, better reproducibility)
```

### 2. Run the Benchmark

```bash
# Run everything
bash main.sh
```

That's it! The system automatically detects your package manager (conda/mamba/micromamba) and environment.

---

## Detailed Setup Options

### Environment Files

Two standardized environment options are provided:

| File | Purpose | When to Use |
|------|---------|------------|
| `environment_minimal.yml` | Flexible versions | Cross-platform compatibility, different systems |
| `environment_complete.yml` | Pinned versions | Exact reproducibility, same platform |

Both create the same environment name: **`benchmarkda`**

### Manual Environment Setup

If you prefer manual setup:

```bash
# Using mamba (recommended - fastest)
mamba env create -f environment_minimal.yml

# Using conda
conda env create -f environment_minimal.yml

# Using micromamba
micromamba env create -f environment_minimal.yml
```

### Package Manager Detection

The system automatically detects and uses:
1. **mamba** (preferred - fastest)
2. **micromamba** (good alternative)
3. **conda** (fallback - slower)

### Adding Kompot from Master

The environment automatically installs kompot from the latest master branch:
```yaml
pip:
  - git+https://github.com/settylab/kompot.git@master
```

---

## Verification

### Test Python Environment

```bash
# Activate environment (auto-detected)
source bin/environment_utils.sh
activate_benchmarkda_environment

# Test core packages
python -c "
import numpy, pandas, scanpy, anndata
import meld, cna, palantir, mellon, kompot
print('✓ All Python DA methods available!')
"
```

### Test R Environment

```bash
Rscript -e "
library(miloR); library(DAseq); library(cydar)
library(SingleCellExperiment); library(scran)
cat('✓ All R DA methods available!\n')
"
```

---

## Troubleshooting

### Environment Not Found

```bash
# Problem: "benchmarkda environment not found"
# Solution: Create the environment
bash setup_environment.sh
```

### Package Manager Not Found

```bash
# Problem: "No conda-compatible package manager found"
# Solution: Install one of these:

# Mamba (recommended)
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh

# Micromamba
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba

# Miniconda
curl -L -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

### Kompot Installation Issues

```bash
# If kompot installation fails, try:
pip install git+https://github.com/settylab/kompot.git@master --force-reinstall

# Or clone and install locally:
git clone https://github.com/settylab/kompot.git
cd kompot
pip install -e .
```

### R Package Issues

```bash
# If R packages are missing:
Rscript -e "
if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')
BiocManager::install(c('miloR', 'SingleCellExperiment', 'scran'))
install.packages(c('argparse', 'tidyverse', 'Seurat', 'igraph'))
"
```

---

## Advanced Usage

### Custom Environment Names

If you need to use a different environment name, you can modify the detection:

```bash
# Set custom environment name
export BENCHMARKDA_ENV_NAME="my_custom_env_name"
bash main.sh
```

### Running Individual Components

```bash
# Generate script for specific dataset/methods
python bin/run_benchmark.py \
    --dataset linear \
    --method_type python \
    --methods kompot mellon \
    --mode_embedding DM \
    --output my_script.sh

# Run the generated script
bash my_script.sh
```

### Platform-Specific Notes

**Linux (Slurm systems)**: Full functionality with parallel job submission
**macOS**: Works, but no Slurm support (sequential execution)
**Windows**: Use WSL2 or Docker

---

## Environment Architecture

The user-agnostic design means:

✅ **No hardcoded paths** - works for any user
✅ **Auto-detection** - finds conda/mamba/micromamba automatically
✅ **Standardized names** - uses `benchmarkda` environment name
✅ **Flexible setup** - supports minimal or complete environments
✅ **Latest kompot** - installs from GitHub master automatically
✅ **Cross-platform** - works on different systems with appropriate environment file

The system replaces user-specific environments like `kompot_v1` with a standardized, reproducible setup that works for anyone.