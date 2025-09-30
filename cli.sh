#!/bin/bash
set -e

# BenchmarkDA: Modular Master Execution Script
# Orchestrates differential abundance benchmarking with granular control

# =============================================================================
# CONFIGURATION AND USAGE
# =============================================================================

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd "$SCRIPT_DIR"

# Source environment utilities
source bin/environment_utils.sh

# Function to load GCC module if needed for R compilation
load_gcc_module_if_needed() {
    # Only load module if module system is available and for R methods
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
                print_info "GLIBC $GLIBC_VERSION detected - loading GCC module for R compilation"
            else
                print_info "GLIBC $GLIBC_VERSION is compatible - no GCC module needed"
            fi
        else
            NEEDS_GCC_MODULE=true
            print_info "Could not detect GLIBC version - loading GCC module as precaution"
        fi

        if [ "$NEEDS_GCC_MODULE" = "true" ]; then
            print_info "Module system detected, attempting to load GCC module for R compilation..."

            # First, completely deactivate any existing mamba/conda environments
            # This ensures module paths will take precedence when we later activate benchmarkda
            print_info "Deactivating any existing conda/mamba environments..."

            # Initialize shell hook if available
            if [ -n "${MAMBA_EXE}" ]; then
                eval "$(${MAMBA_EXE} shell hook --shell bash)" 2>/dev/null
                while [ ! -z "$CONDA_PREFIX" ]; do
                    mamba deactivate 2>/dev/null || break
                done
            elif [ -n "${CONDA_EXE}" ]; then
                eval "$(${CONDA_EXE} shell hook --shell bash)" 2>/dev/null
                while [ ! -z "$CONDA_PREFIX" ]; do
                    conda deactivate 2>/dev/null || break
                done
            elif command -v micromamba &> /dev/null; then
                eval "$(micromamba shell hook --shell bash)" 2>/dev/null
                while [ ! -z "$CONDA_PREFIX" ]; do
                    micromamba deactivate 2>/dev/null || break
                done
            fi

            print_info "All environments deactivated before module load"

            # Load GCC module (will be at front of PATH)
            for gcc_version in "13.3.0" "13.2.0" "12.3.0" "12.2.0" "11.3.0" "11.2.0"; do
                if module avail GCC/$gcc_version 2>&1 | grep -q "GCC/$gcc_version"; then
                    print_info "Loading GCC/$gcc_version module for R compilation support"
                    module load GCC/$gcc_version
                    return 0
                fi
            done
            print_warning "No compatible GCC modules found"
        fi
    fi
}

# Color codes for better output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m'

# Configuration - hardcoded for shell compatibility
SYNTHETIC_DATASETS=(linear branch cluster)
REAL_DATASETS=(covid19-pbmc pancreas bcr-xl levine32)

print_header() {
    echo -e "${BOLD}${BLUE}========================================${NC}"
    echo -e "${BOLD}${BLUE}$1${NC}"
    echo -e "${BOLD}${BLUE}========================================${NC}"
}

print_step() {
    echo -e "${GREEN}[STEP]${NC} $1"
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

usage() {
    cat << EOF
BenchmarkDA: Differential Abundance Benchmarking Pipeline

USAGE: ./cli.sh [OPTIONS] [STEPS...]

OPTIONS:
    -h, --help              Show this help message
    -d, --datasets LIST     Datasets: linear,branch,cluster,covid19-pbmc,pancreas,bcr-xl,levine32
    -m, --methods LIST      Methods: python,r,all (default: all)
    -s, --skip-missing      Skip missing datasets without prompting
    --slurm                Submit benchmarks as SLURM array jobs
    --dry-run              Show commands without executing

STEPS:
    setup                   Setup environment and directories
    preprocess             Create PCA embeddings and DM from PCA+batch
    labels                 Generate synthetic labels
    benchmark              Run DA method benchmarks
    status                 Show completion status
    all                    Run all steps (default)

STATUS COMMANDS:
    ./cli.sh status                              All datasets
    ./cli.sh --datasets linear,branch status     Specific datasets

STATUS OUTPUT:
    [Complete] = Finished   [Partial] = In progress   [Missing] = Not started

PYTHON METHODS:
    mellon, mellon_noSync, mellon_corr, mellon_pca
    meld, meld_default, meld_pca
    kompot, kompot_pca

R METHODS:
    milo, daseq, cydar, louvain

EXAMPLES:
    ./cli.sh                                     Complete pipeline (local)
    ./cli.sh benchmark --slurm                   Submit ALL benchmarks to SLURM
    ./cli.sh status                              Check progress
    ./cli.sh --datasets linear preprocess        Preprocess one dataset
    ./cli.sh --methods python benchmark          Python methods only
    ./cli.sh --methods r benchmark --slurm       Submit only R methods to SLURM
    ./cli.sh --dry-run                           Show commands without execution
EOF
}

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

# Default values
DATASETS=""
METHODS="all"
SKIP_MISSING=false
DRY_RUN=false
STEPS=()

while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            usage
            exit 0
            ;;
        -d|--datasets)
            DATASETS="$2"
            shift 2
            ;;
        -m|--methods)
            METHODS="$2"
            shift 2
            ;;
        -s|--skip-missing)
            SKIP_MISSING=true
            shift
            ;;
        --slurm)
            USE_SLURM=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        setup|preprocess|labels|benchmark|status|all)
            STEPS+=("$1")
            shift
            ;;
        *)
            echo "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Default to all steps if none specified
if [ ${#STEPS[@]} -eq 0 ]; then
    STEPS=("all")
fi

# Parse datasets
if [ -z "$DATASETS" ]; then
    SELECTED_DATASETS=(${SYNTHETIC_DATASETS[@]} ${REAL_DATASETS[@]})
else
    IFS=',' read -ra SELECTED_DATASETS <<< "$DATASETS"
fi

# Parse methods
if [ "$METHODS" = "all" ]; then
    SELECTED_METHODS=("python" "r")
else
    IFS=',' read -ra SELECTED_METHODS <<< "$METHODS"
fi

# =============================================================================
# CORE FUNCTIONS
# =============================================================================

setup_environment() {
    print_step "Setting up computational environment"

    # Load environment utilities
    if [[ -f "bin/environment_utils.sh" ]]; then
        source "bin/environment_utils.sh"
        print_info "Loaded environment detection utilities"
    else
        print_error "Environment utilities not found!"
        exit 1
    fi

    # Load GCC module first if needed (before mamba activation)
    load_gcc_module_if_needed

    # Check if benchmarkda environment exists
    ENV_NAME=$(detect_benchmarkda_environment)
    if [[ -z "$ENV_NAME" ]]; then
        print_warning "BenchmarkDA environment not found!"
        if [ "$SKIP_MISSING" = false ]; then
            echo "Would you like to create it now? [y/N]"
            read -r response
            if [[ "$response" =~ ^[Yy]$ ]]; then
                print_step "Creating benchmarkda environment"
                bash setup_environment.sh --minimal
                ENV_NAME=$(detect_benchmarkda_environment)
            fi
        fi

        if [[ -z "$ENV_NAME" ]]; then
            print_error "Cannot proceed without benchmarkda environment"
            print_error "Run: bash setup_environment.sh"
            exit 1
        fi
    fi

    # Activate environment using shell hook (sets up environment for inheritance)
    print_info "Activating environment: $ENV_NAME"
    if ! activate_benchmarkda_environment; then
        print_error "Failed to activate environment"
        exit 1
    fi

    # Create directory structure
    print_info "Creating directory structure"
    mkdir -p {data/{synthetic,real/{bcr-xl,covid19-pbmc,levine32,pancreas}},SlurmLog,benchmark_scripts}

    print_success "Environment setup completed"
}

get_n_dm_for_dataset() {
    local dataset="$1"

    # Try Python config first, fallback to hardcoded values
    if python -c "
import sys
sys.path.append('.')
try:
    from config.dataset_config import get_n_dm_for_dataset
    print(get_n_dm_for_dataset('$dataset'))
except:
    pass
" 2>/dev/null | grep -q '^[0-9]*$'; then
        python -c "
import sys
sys.path.append('.')
from config.dataset_config import get_n_dm_for_dataset
print(get_n_dm_for_dataset('$dataset'))
"
    else
        # Fallback hardcoded values
        case "$dataset" in
            linear|branch|cluster)
                echo "10"
                ;;
            covid19-pbmc|pancreas)
                echo "30"
                ;;
            bcr-xl|levine32)
                echo "5"
                ;;
            *)
                echo "10"  # default
                ;;
        esac
    fi
}

check_datasets() {
    print_step "Checking dataset availability"

    local missing_datasets=()
    for dataset in "${SELECTED_DATASETS[@]}"; do
        if [[ " ${SYNTHETIC_DATASETS[*]} " =~ " ${dataset} " ]]; then
            # Synthetic dataset - should exist or be generated
            continue
        elif [[ " ${REAL_DATASETS[*]} " =~ " ${dataset} " ]]; then
            # Real dataset - check if file exists
            if [[ ! -f "data/real/$dataset/$dataset.h5ad" ]]; then
                missing_datasets+=("$dataset")
            fi
        else
            print_warning "Unknown dataset: $dataset"
        fi
    done

    if [[ ${#missing_datasets[@]} -gt 0 ]]; then
        print_warning "Missing real datasets: ${missing_datasets[*]}"
        print_info "Download from: https://drive.google.com/drive/folders/15wWFD5FMe0VdzN1pUnaUUpQ17OXkeebH"

        if [ "$SKIP_MISSING" = false ]; then
            echo "Continue without real datasets? [Y/n]"
            read -r response
            if [[ "$response" =~ ^[Nn]$ ]]; then
                exit 1
            fi
        fi

        # Remove missing datasets from selection
        for missing in "${missing_datasets[@]}"; do
            SELECTED_DATASETS=(${SELECTED_DATASETS[@]/$missing})
        done
    fi

    print_info "Processing datasets: ${SELECTED_DATASETS[*]}"
}

preprocess_datasets() {
    print_step "Preprocessing datasets"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        local n_dm=$(get_n_dm_for_dataset "$dataset")

        print_info "Processing $dataset (DM components: $n_dm)"

        # First create PCA embeddings
        if [ "$DRY_RUN" = true ]; then
            echo "[DRY RUN] bash bin/dataset_preprocessing.sh \"$dataset\" X_pca 0 PCA"
        else
            bash bin/dataset_preprocessing.sh "$dataset" X_pca 0 PCA || {
                print_warning "PCA preprocessing failed for $dataset"
                continue
            }
        fi

        # Then create DM from PCA (for methods that need it)
        if [[ $n_dm -gt 0 ]]; then
            if [ "$DRY_RUN" = true ]; then
                echo "[DRY RUN] bash bin/dataset_preprocessing.sh \"$dataset\" X_pca \"$n_dm\" DM"
            else
                bash bin/dataset_preprocessing.sh "$dataset" X_pca "$n_dm" DM || {
                    print_warning "DM preprocessing failed for $dataset"
                }
            fi
        fi
    done

    print_success "Dataset preprocessing completed"
}

generate_labels() {
    print_step "Generating synthetic condition labels"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        print_info "Generating labels for $dataset"

        if [ "$DRY_RUN" = true ]; then
            echo "[DRY RUN] bash bin/modified_benchmarkda_dm_all.sh \"$dataset\" labels No PCA 0"
        else
            # Generate labels once per dataset
            bash bin/modified_benchmarkda_dm_all.sh "$dataset" labels No PCA 0 || {
                print_warning "Label generation failed for $dataset"
            }
        fi
    done

    print_success "Label generation completed"
}

run_benchmarks() {
    print_step "Running differential abundance method benchmarks"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        local n_dm=$(get_n_dm_for_dataset "$dataset")

        for method_type in "${SELECTED_METHODS[@]}"; do
            print_info "Running $method_type methods on $dataset with embeddings: ${SELECTED_EMBEDDINGS[*]}"

            # Prepare SLURM flag if needed
            local slurm_flag=""
            if [ "$USE_SLURM" = true ]; then
                slurm_flag="--slurm"
            fi

            if [ "$DRY_RUN" = true ]; then
                echo "[DRY RUN] python bin/direct_benchmark.py --dataset \"$dataset\" --method_type \"$method_type\" --n_dm \"$n_dm\" $slurm_flag"
            else
                if [ "$USE_SLURM" = true ]; then
                    print_info "Submitting $method_type methods for $dataset to SLURM"
                else
                    print_info "Executing benchmarks"
                fi

                # Load GCC module if needed for R methods before mamba activation
                if [ "$method_type" = "r" ]; then
                    load_gcc_module_if_needed
                fi

                # Use proper environment execution command
                local run_cmd=$(get_environment_run_command)
                if [ -n "$run_cmd" ]; then
                    $run_cmd python bin/direct_benchmark.py \
                        --dataset "$dataset" \
                        --method_type "$method_type" \
                        --n_dm "$n_dm" \
                        $slurm_flag || {
                        print_warning "Direct benchmark execution failed for $dataset $method_type"
                    }
                else
                    print_warning "No suitable environment execution command found, trying direct execution"
                    python bin/direct_benchmark.py \
                        --dataset "$dataset" \
                        --method_type "$method_type" \
                        --n_dm "$n_dm" \
                        $slurm_flag || {
                        print_warning "Direct benchmark execution failed for $dataset $method_type"
                    }
                fi

                # If using SLURM, submit the generated script
                if [ "$USE_SLURM" = true ] && [ "$DRY_RUN" != true ]; then
                    local script_path="benchmark_scripts/slurm_benchmark_${dataset}_${method_type}.sh"
                    if [ -f "$script_path" ]; then
                        local job_id=$(sbatch --export=ALL "$script_path" | awk '{print $NF}')
                        print_success "Submitted SLURM job $job_id for $dataset $method_type"
                    else
                        print_warning "SLURM script not found: $script_path"
                    fi
                fi
            fi
        done
    done

    print_success "Benchmark execution completed"
}

show_status_report() {
    print_step "Generating status report"

    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] python bin/status_report.py --detailed --commands --save"
    else
        # Use MAMBA_EXE if available for reliable environment execution
        if [ -n "${MAMBA_EXE}" ]; then
            if [ ${#SELECTED_DATASETS[@]} -gt 0 ]; then
                # Report on specific datasets
                for dataset in "${SELECTED_DATASETS[@]}"; do
                    print_info "Status for dataset: $dataset"
                    ${MAMBA_EXE} run -n benchmarkda python bin/status_report.py --dataset "$dataset" --detailed
                    echo
                done
            else
                # Report on all datasets
                ${MAMBA_EXE} run -n benchmarkda python bin/status_report.py --detailed --commands --save
            fi
        else
            print_warning "MAMBA_EXE not available, trying direct execution"
            if [ ${#SELECTED_DATASETS[@]} -gt 0 ]; then
                for dataset in "${SELECTED_DATASETS[@]}"; do
                    print_info "Status for dataset: $dataset"
                    python bin/status_report.py --dataset "$dataset" --detailed
                    echo
                done
            else
                python bin/status_report.py --detailed --commands --save
            fi
        fi
    fi

    print_success "Status report completed"
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

main() {
    print_header "BenchmarkDA: Modular Differential Abundance Benchmarking"

    if [ "$DRY_RUN" = true ]; then
        print_info "DRY RUN MODE - No actual execution"
    fi

    print_info "Configuration:"
    print_info "  Datasets: ${SELECTED_DATASETS[*]}"
    print_info "  Methods: ${SELECTED_METHODS[*]}"
    print_info "  Embeddings: PCA-based"
    print_info "  Steps: ${STEPS[*]}"

    # Check if we should run all steps
    if [[ " ${STEPS[*]} " =~ " all " ]]; then
        STEPS=("setup" "preprocess" "labels" "benchmark")
    fi

    # Execute steps
    for step in "${STEPS[@]}"; do
        case $step in
            setup)
                setup_environment
                check_datasets
                ;;
            preprocess)
                preprocess_datasets
                ;;
            labels)
                generate_labels
                ;;
            benchmark)
                run_benchmarks
                ;;
            status)
                show_status_report
                ;;
            *)
                print_error "Unknown step: $step"
                exit 1
                ;;
        esac
    done

    # Summary
    print_header "Execution Summary"
    echo ""
    print_success "COMPLETED STEPS: ${STEPS[*]}"
    print_info "Results location: benchmark/"
    print_info "Generated scripts: benchmark_scripts/"
    print_info "Slurm logs: SlurmLog/"
    echo ""
    print_step "INDIVIDUAL EXECUTION EXAMPLES"
    echo "# Run specific method on dataset:"
    echo "python run_da_method.py kompot linear --embedding dm"
    echo ""
    echo "# Run specific components:"
    echo "./cli.sh --datasets linear --methods python preprocess benchmark"

    print_header "All Done!"
}

main "$@"