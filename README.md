# BenchmarkDA: Differential Abundance Testing Framework

A comprehensive, user-agnostic framework for benchmarking differential abundance (DA) methods in single-cell data with consistent batch-corrected embeddings.

## 🚀 Quick Start

```bash
# 1. Test your setup
python test_setup.py

# 2. Create environment (if needed)
bash setup_environment.sh

# 3. Run complete benchmark
bash main.sh

# 4. Run individual components
./cli.sh --datasets linear,branch preprocess
./cli.sh --methods python --embeddings dm benchmark

# 5. Run individual methods
python run_da_method.py kompot linear --embedding dm
```

## 📋 What's Included

### Python Methods
- **kompot** (from master): Latest from `github.com/settylab/kompot`
- **mellon**: Multiple variants (standard, PCA, normalized, corrected)
- **meld**: Multiple variants (standard, PCA, default)
- **cna**: Multiple variants (standard, PCA)

### R Methods
- **milo**: Neighborhood-based DA testing
- **daseq**: DA region detection
- **cydar**: Hypersphere-based DA testing
- **louvain**: Clustering-based DA testing

### Datasets
- **Synthetic**: Linear, branch, cluster topologies
- **Real**: BCR-XL, COVID-19 PBMC, Levine32, pancreas (download separately)

## 🏗️ Architecture Highlights

### ✅ **Fixed Critical Issues**
- **Batch Effect Consistency**: All methods use same `DM_EigenVectors_batch` and `X_pca_batch`
- **User-Agnostic**: Works with any conda/mamba/micromamba setup
- **Modular Design**: Run individual methods or complete pipeline
- **Clean Structure**: Organized, readable, maintainable code

### ✅ **Improved Usability**
- **One-Command Setup**: `bash setup_environment.sh`
- **Comprehensive Testing**: `python test_setup.py`
- **Clear Main Script**: Readable `main.sh` with colored output and error handling
- **Individual Execution**: `python run_da_method.py` for single methods

## 📂 Directory Structure

```
benchmarkDA_private/
├── main.sh                          # 🎯 Master workflow (clean & readable)
├── setup_environment.sh             # 🔧 User-agnostic environment setup
├── test_setup.py                    # ✅ Comprehensive testing
├── run_da_method.py                 # 🎮 Individual method execution
│
├── config/
│   ├── method_config.py             # Method configurations
│   └── dataset_config.py            # Dataset parameters
│
├── python_method/
│   ├── shared_embedding_utils.py    # 🔑 Critical: Consistent embeddings
│   ├── data_loader.py               # Unified data loading
│   ├── Mellon_bm.py, meld_bm.py     # Method implementations
│   ├── CNA_bm.py, kompot_bm.py      # (all updated for consistency)
│   └── run*.py                      # Method-specific runners
│
├── scripts/
│   └── run_DA.r                     # 🔄 Unified R methods
│
├── bin/
│   ├── environment_utils.sh         # Environment detection
│   └── run_benchmark.py             # Script generator
│
├── environment_minimal.yml          # Flexible versions
├── environment_complete.yml         # Pinned versions
└── legacy/                          # Old files moved here
```

## 🎯 Usage Examples

### Complete Pipeline
```bash
# Full benchmark with all methods and datasets
bash main.sh
```

### Modular Execution
```bash
# List available methods and datasets
python run_da_method.py --list

# Run specific pipeline steps
./cli.sh --datasets linear,branch preprocess
./cli.sh --methods python --embeddings dm benchmark
./cli.sh --skip-missing --dry-run all

# Run individual methods
python run_da_method.py kompot linear --embedding dm
python run_da_method.py mellon branch --embedding pca --population M2
python run_da_method.py meld cluster --embedding dm --enrichment 3.0 --seed 42
```

### Generated Scripts (for advanced users)
```bash
# Generate script for specific dataset/methods
python bin/run_benchmark.py \
    --dataset linear \
    --method_type python \
    --methods kompot mellon \
    --mode_embedding DM

# Generated scripts are saved in benchmark_scripts/
```

## 🔧 Environment Management

### Automatic Setup
```bash
bash setup_environment.sh          # Choose minimal or complete
```

### Manual Setup
```bash
# Option 1: Minimal (flexible versions)
mamba env create -f environment_minimal.yml

# Option 2: Complete (pinned versions)
mamba env create -f environment_complete.yml
```

### Environment Detection
The system automatically detects:
- **Package managers**: mamba → micromamba → conda
- **Environment names**: benchmarkda, diffabundance, kompot_v1, etc.
- **Missing packages**: Provides clear error messages

## ✅ Testing & Validation

```bash
# Comprehensive test suite
python test_setup.py

# Test specific components
python test_setup.py  # Will test:
# ✓ File structure
# ✓ Method configurations
# ✓ Environment detection
# ✓ Python packages (including kompot from master)
# ✓ R packages
# ✓ Basic functionality
```

## 🚨 What Was Fixed

### Old Issues ❌
- Hardcoded user-specific environment names (`kompot_v1`)
- Inconsistent diffusion map computation across methods
- Hard-to-read monolithic `main.sh` script
- No modular execution options
- Cluttered repository with old files
- Complex script generation system

### New Solutions ✅
- **User-agnostic environment** (`benchmarkda`) with auto-detection
- **Consistent embeddings** via `shared_embedding_utils.py`
- **Clean, readable main.sh** with colored output and error handling
- **Modular execution** via `run_da_method.py`
- **Organized structure** with `legacy/` folder for old files
- **Simple testing** with `test_setup.py`

## 📊 Output Structure

```
benchmark/                          # Results organized by embedding
├── dm/                             # Diffusion map results
│   ├── synthetic/{linear,branch,cluster}/
│   └── real/{covid19-pbmc,pancreas,...}/
└── pca/                           # PCA results (same structure)

results/individual/                 # Individual method results
└── {method}_{dataset}_{embedding}/

benchmark_scripts/                  # Generated scripts (for inspection)
└── {dataset}_{embedding}_{type}.sh
```

## 🎯 Why This Is Better

| Aspect | Before | After |
|--------|---------|--------|
| **Environment** | User-specific (`kompot_v1`) | User-agnostic (`benchmarkda`) |
| **Embedding Consistency** | ❌ Inconsistent DM computation | ✅ Shared embedding utilities |
| **Modularity** | ❌ Monolithic scripts only | ✅ Individual method execution |
| **Readability** | ❌ Hard-to-read main.sh | ✅ Clean, colored, modular |
| **Testing** | ❌ No validation system | ✅ Comprehensive test suite |
| **Organization** | ❌ Cluttered with old files | ✅ Clean structure + legacy/ |
| **Setup** | ❌ Manual, error-prone | ✅ One-command setup + validation |

## 🔄 Migration from Old Setup

If you had the old `kompot_v1` setup:

```bash
# Your old setup still works! But to get the improvements:
git pull                          # Get latest changes
python test_setup.py             # Test your current setup
bash setup_environment.sh        # Create standardized environment
bash main.sh                     # Run with new architecture
```

---

**Ready to benchmark differential abundance methods reliably!** 🎉