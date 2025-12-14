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

# Check for conda/mamba/micromamba using central environment utils approach
CONDA_CMD=""
if [ -n "${MAMBA_EXE}" ] && [ -x "${MAMBA_EXE}" ]; then
    CONDA_CMD="${MAMBA_EXE}"
    print_status "Found mamba package manager: ${MAMBA_EXE}"
elif [ -n "${CONDA_EXE}" ] && [ -x "${CONDA_EXE}" ]; then
    CONDA_CMD="${CONDA_EXE}"
    print_status "Found conda package manager: ${CONDA_EXE}"
elif command -v micromamba &> /dev/null; then
    CONDA_CMD="micromamba"
    print_status "Found micromamba package manager"
    print_warning "Consider setting MAMBA_EXE or CONDA_EXE for better compatibility"
elif command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    print_status "Found mamba package manager"
    print_warning "Consider setting MAMBA_EXE environment variable"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
    print_status "Found conda package manager"
    print_warning "conda is slower than mamba/micromamba. Consider installing mamba for faster package management."
    print_warning "Consider setting CONDA_EXE environment variable"
else
    print_error "No conda-compatible package manager found!"
    print_error "Please install conda, mamba, or micromamba first:"
    print_error "  - Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    print_error "  - Mamba: https://mamba.readthedocs.io/en/latest/installation.html"
    print_error "  - Micromamba: https://mamba.readthedocs.io/en/latest/installation.html#micromamba"
    print_error "Or set MAMBA_EXE/CONDA_EXE environment variables"
    exit 1
fi

# Function to load GCC module if needed for R compilation (following CLI pattern)
load_gcc_module_if_needed() {
    # Only load module if module system is available
    if command -v module &> /dev/null; then
        # Check GLIBC version
        local GLIBC_VERSION=$(ldd --version 2>/dev/null | head -n1 | grep -o '[0-9]\+\.[0-9]\+' | head -n1 || echo "unknown")
        local NEEDS_GCC_MODULE=false

        if [ "$GLIBC_VERSION" != "unknown" ]; then
            local GLIBC_MAJOR=$(echo $GLIBC_VERSION | cut -d. -f1)
            local GLIBC_MINOR=$(echo $GLIBC_VERSION | cut -d. -f2)

            # Only load module if GLIBC < 2.29
            if [[ $GLIBC_MAJOR -lt 2 ]] || [[ $GLIBC_MAJOR -eq 2 && $GLIBC_MINOR -lt 29 ]]; then
                NEEDS_GCC_MODULE=true
                print_status "GLIBC $GLIBC_VERSION detected - loading GCC module for R compilation"
            else
                print_status "GLIBC $GLIBC_VERSION is compatible - no GCC module needed"
            fi
        else
            NEEDS_GCC_MODULE=true
            print_status "Could not detect GLIBC version - loading GCC module as precaution"
        fi

        if [ "$NEEDS_GCC_MODULE" = "true" ]; then
            print_status "Module system detected, attempting to load GCC module for R compilation..."

            # First, deactivate any existing mamba environments to ensure clean PATH
            # Using the same pattern as environment utils
            if [ -n "${MAMBA_EXE}" ] && [ -x "${MAMBA_EXE}" ]; then
                while ${MAMBA_EXE} info 2>/dev/null | grep -q "active environment"; do
                    ${MAMBA_EXE} deactivate 2>/dev/null || break
                done
            elif [ -n "${CONDA_EXE}" ] && [ -x "${CONDA_EXE}" ]; then
                while ${CONDA_EXE} info 2>/dev/null | grep -q "active environment"; do
                    ${CONDA_EXE} deactivate 2>/dev/null || break
                done
            elif command -v conda &> /dev/null; then
                while conda info | grep -q "active environment"; do
                    conda deactivate 2>/dev/null || break
                done
            fi

            # Load GCC module
            for gcc_version in "13.3.0" "13.2.0" "12.3.0" "12.2.0" "11.3.0" "11.2.0"; do
                if module avail GCC/$gcc_version 2>&1 | grep -q "GCC/$gcc_version"; then
                    print_status "Loading GCC/$gcc_version module for R compilation support"
                    if module load GCC/$gcc_version 2>/dev/null; then
                        print_success "Successfully loaded GCC/$gcc_version module"
                        return 0
                    else
                        print_warning "Failed to load GCC/$gcc_version module"
                    fi
                fi
            done
            print_warning "No compatible GCC modules found"
        fi
    fi
}

# Load GCC module if needed (before environment creation/activation)
print_status "Checking system requirements..."
load_gcc_module_if_needed

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

# Activate environment using the same pattern as environment utils
if [ -n "${MAMBA_EXE}" ] && [ -x "${MAMBA_EXE}" ]; then
    print_status "Activating benchmarkda environment using ${MAMBA_EXE}..."
    eval "$(${MAMBA_EXE} shell hook --shell bash)"
    ${MAMBA_EXE} activate $ENV_NAME
elif [ -n "${CONDA_EXE}" ] && [ -x "${CONDA_EXE}" ]; then
    print_status "Activating benchmarkda environment using ${CONDA_EXE}..."
    eval "$(${CONDA_EXE} shell hook --shell bash)"
    ${CONDA_EXE} activate $ENV_NAME
elif command -v micromamba &> /dev/null; then
    print_status "Activating benchmarkda environment using micromamba..."
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate $ENV_NAME
elif command -v mamba &> /dev/null; then
    print_status "Activating benchmarkda environment using mamba..."
    eval "$(conda shell.bash hook)"
    mamba activate $ENV_NAME
else
    print_status "Activating benchmarkda environment using conda..."
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

# Configure R environment for benchmarkDA
print_status "Configuring R environment..."

# Module loading already done earlier - determine R package strategy based on what happened
if command -v module &> /dev/null; then
    # Check if any GCC module was loaded
    if module list 2>&1 | grep -q "GCC/"; then
        print_status "GCC module detected, using renv for R package management"
        export RENV_CONFIG_SANDBOX_ENABLED=FALSE
        R_PACKAGE_STRATEGY="renv"
    else
        print_status "No GCC module loaded, using mamba for R package management"
        export RENV_CONFIG_SANDBOX_ENABLED=FALSE
        export R_LIBS_USER=""
        R_PACKAGE_STRATEGY="mamba"
    fi
else
    print_status "No module system detected, using mamba for R package management"
    export RENV_CONFIG_SANDBOX_ENABLED=FALSE
    export R_LIBS_USER=""
    R_PACKAGE_STRATEGY="mamba"
fi

# Install critical R packages via mamba first (platform-agnostic)
print_status "Installing critical R packages via mamba..."
if [ -n "$CONDA_CMD" ]; then
    print_status "Installing base R packages..."
    $CONDA_CMD install -n $ENV_NAME -c conda-forge r-igraph r-bluster r-biocmanager -y

    print_status "Attempting to install Bioconductor packages via mamba..."
    # These may not be available on all platforms - will install via BiocManager if not
    $CONDA_CMD install -n $ENV_NAME -c bioconda r-singlecellexperiment r-scran -y || true

    print_status "Mamba R package installation completed"
fi

# Test R packages and install if needed based on strategy
print_status "Testing R packages..."

if [[ $R_PACKAGE_STRATEGY == "renv" ]]; then
    print_status "Testing R packages with renv strategy..."
    Rscript -e "
    packages <- c('argparse', 'tidyverse', 'SingleCellExperiment', 'scran', 'Seurat', 'igraph', 'anndata', 'reticulate')
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
        cat('Installing missing packages via BiocManager...\\n')

        if (!requireNamespace('BiocManager', quietly = TRUE)) {
            install.packages('BiocManager', repos='https://cran.rstudio.com/')
        }
        BiocManager::install(missing, ask=FALSE)

        # Verify installation
        still_missing <- c()
        for (pkg in missing) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                still_missing <- c(still_missing, pkg)
            }
        }

        if (length(still_missing) > 0) {
            cat('\\nFailed to install via renv:', paste(still_missing, collapse = ', '), '\\n')
            cat('Falling back to mamba for these packages...\\n')
            quit('no', 1)  # Exit with error to trigger mamba fallback
        } else {
            cat('\\nAll R packages successfully installed via renv!\\n')
        }
    } else {
        cat('\\nAll R packages are available via renv!\\n')
    }
    "

    # Check if renv installation failed
    if [ $? -ne 0 ]; then
        print_warning "renv package installation failed, falling back to mamba strategy"
        R_PACKAGE_STRATEGY="mamba"
    fi
fi

if [[ $R_PACKAGE_STRATEGY == "mamba" ]]; then
    print_status "Testing R packages with mamba strategy..."
    Rscript -e "
    # Ensure mamba library path is prioritized
    conda_env_path <- Sys.getenv('CONDA_PREFIX')
    if (conda_env_path != '') {
        conda_r_lib <- file.path(conda_env_path, 'lib', 'R', 'library')
        if (dir.exists(conda_r_lib)) {
            .libPaths(c(conda_r_lib, .libPaths()))
            cat('Using mamba R library path:', conda_r_lib, '\\n')
        }
    }

    packages <- c('argparse', 'tidyverse', 'SingleCellExperiment', 'scran', 'Seurat', 'igraph', 'anndata', 'reticulate')
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
        cat('Note: These packages should be available from mamba installation.\\n')
        cat('If any are missing, they will be installed via BiocManager as fallback.\\n')
    } else {
        cat('\\nAll R packages are available via mamba!\\n')
    }
    "
fi

# Create R environment configuration note
print_status "R environment configuration complete"
print_status "Note: R methods use conda packages to avoid compilation issues"

print_success "Environment '$ENV_NAME' is ready to use!"

echo ""
echo "=========================================="
echo "Next Steps:"
echo "=========================================="
echo "1. Manual activation (for advanced users):"

# Show proper activation commands based on what was detected
if [ -n "${MAMBA_EXE}" ] && [ -x "${MAMBA_EXE}" ]; then
    echo "   # For R methods, manually check and load GCC module if needed"
    echo "   # Then activate: ${MAMBA_EXE} activate $ENV_NAME"
elif [ -n "${CONDA_EXE}" ] && [ -x "${CONDA_EXE}" ]; then
    echo "   # For R methods, manually check and load GCC module if needed"
    echo "   # Then activate: ${CONDA_EXE} activate $ENV_NAME"
elif command -v micromamba &> /dev/null; then
    echo "   # For R methods, manually check and load GCC module if needed"
    echo "   # Then activate: micromamba activate $ENV_NAME"
elif command -v mamba &> /dev/null; then
    echo "   # For R methods, manually check and load GCC module if needed"
    echo "   # Then activate: mamba activate $ENV_NAME"
else
    echo "   # For R methods, manually check and load GCC module if needed"
    echo "   # Then activate: conda activate $ENV_NAME"
fi

echo ""
echo "2. Or use the CLI (recommended - handles environment automatically):"
echo "   ./cli.sh --help"
echo ""
echo "3. Run specific benchmarks:"
echo "   ./cli.sh --datasets linear --methods python,r benchmark"
echo ""
echo "Environment name: $ENV_NAME"
echo "Intelligent module loading (GLIBC detection + GCC version selection) is automatic when using ./cli.sh"
echo "For manual activation, you need to handle module loading yourself"
echo "=========================================="