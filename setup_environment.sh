#!/bin/bash

# BenchmarkDA Environment Setup Script
# Creates a user-agnostic environment for running differential abundance benchmarks

set -e

echo "=========================================="
echo "BenchmarkDA Environment Setup"
echo "=========================================="

# Color codes for better output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check for conda/mamba/micromamba
CONDA_CMD=""
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    print_status "Found mamba package manager"
elif command -v micromamba &> /dev/null; then
    CONDA_CMD="micromamba"
    print_status "Found micromamba package manager"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    print_status "Found conda package manager"
    print_warning "conda is slower than mamba/micromamba. Consider installing mamba for faster package management."
else
    print_error "No conda-compatible package manager found!"
    print_error "Please install conda, mamba, or micromamba first:"
    print_error "  - Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    print_error "  - Mamba: https://mamba.readthedocs.io/en/latest/installation.html"
    print_error "  - Micromamba: https://mamba.readthedocs.io/en/latest/installation.html#micromamba"
    exit 1
fi

# Choose environment file
ENV_FILE=""
if [ "$1" = "--minimal" ]; then
    ENV_FILE="environment_minimal.yml"
    print_status "Using minimal environment (flexible versions)"
elif [ "$1" = "--complete" ]; then
    ENV_FILE="environment_complete.yml"
    print_status "Using complete environment (pinned versions)"
else
    print_status "Choose environment type:"
    echo "  1) Minimal (flexible versions, better compatibility)"
    echo "  2) Complete (pinned versions, better reproducibility)"
    read -p "Enter choice [1-2]: " choice

    case $choice in
        1)
            ENV_FILE="environment_minimal.yml"
            print_status "Using minimal environment"
            ;;
        2)
            ENV_FILE="environment_complete.yml"
            print_status "Using complete environment"
            ;;
        *)
            print_warning "Invalid choice, defaulting to minimal environment"
            ENV_FILE="environment_minimal.yml"
            ;;
    esac
fi

# Check if environment file exists
if [ ! -f "$ENV_FILE" ]; then
    print_error "Environment file $ENV_FILE not found!"
    exit 1
fi

# Check if environment already exists
ENV_NAME="benchmarkda"
if $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
    print_warning "Environment '$ENV_NAME' already exists!"
    read -p "Do you want to update it? [y/N]: " update_env
    if [[ $update_env =~ ^[Yy]$ ]]; then
        print_status "Updating existing environment..."
        $CONDA_CMD env update -f "$ENV_FILE"
    else
        print_status "Skipping environment creation"
    fi
else
    print_status "Creating new environment '$ENV_NAME'..."
    $CONDA_CMD env create -f "$ENV_FILE"
fi

print_success "Environment setup completed!"

# Test the environment
print_status "Testing environment..."
if command -v micromamba &> /dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate $ENV_NAME
elif command -v mamba &> /dev/null; then
    eval "$(conda shell.bash hook)"
    mamba activate $ENV_NAME
else
    eval "$(conda shell.bash hook)"
    conda activate $ENV_NAME
fi

# Test Python packages
print_status "Testing Python packages..."
python -c "
import sys
packages = ['numpy', 'pandas', 'scanpy', 'anndata', 'meld', 'cna', 'palantir', 'mellon']
missing = []
for pkg in packages:
    try:
        __import__(pkg)
        print(f'✓ {pkg}')
    except ImportError:
        missing.append(pkg)
        print(f'✗ {pkg}')

# Test kompot specifically
try:
    import kompot
    print('✓ kompot (from GitHub master)')
except ImportError:
    missing.append('kompot')
    print('✗ kompot')

if missing:
    print(f'\\nMissing packages: {missing}')
    print('You may need to install them manually or check your environment.')
    sys.exit(1)
else:
    print('\\nAll Python packages are available!')
"

# Test R packages
print_status "Testing R packages..."
Rscript -e "
packages <- c('argparse', 'tidyverse', 'SingleCellExperiment', 'scran', 'Seurat', 'igraph')
missing <- c()
for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        missing <- c(missing, pkg)
        cat(paste('✗', pkg, '\\n'))
    } else {
        cat(paste('✓', pkg, '\\n'))
    }
}

if (length(missing) > 0) {
    cat('\\nMissing R packages:', paste(missing, collapse = ', '), '\\n')
    cat('Installing missing R packages...\\n')

    # Install Bioconductor packages
    if (!requireNamespace('BiocManager', quietly = TRUE)) {
        install.packages('BiocManager', repos = 'https://cloud.r-project.org/')
    }

    BiocManager::install(missing[missing %in% c('SingleCellExperiment', 'scran')])

    # Install CRAN packages
    cran_packages <- missing[!missing %in% c('SingleCellExperiment', 'scran')]
    if (length(cran_packages) > 0) {
        install.packages(cran_packages, repos = 'https://cloud.r-project.org/')
    }
} else {
    cat('\\nAll R packages are available!\\n')
}
"

print_success "Environment '$ENV_NAME' is ready to use!"

echo ""
echo "=========================================="
echo "Next Steps:"
echo "=========================================="
echo "1. Activate the environment:"
if command -v micromamba &> /dev/null; then
    echo "   micromamba activate $ENV_NAME"
elif command -v mamba &> /dev/null; then
    echo "   mamba activate $ENV_NAME"
else
    echo "   conda activate $ENV_NAME"
fi
echo ""
echo "2. Run the benchmark:"
echo "   bash main.sh"
echo ""
echo "3. Or run individual methods:"
echo "   python bin/run_benchmark.py --help"
echo ""
echo "Environment name: $ENV_NAME"
echo "This name is used throughout the benchmarking scripts."
echo "=========================================="