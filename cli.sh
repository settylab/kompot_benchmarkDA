#!/bin/bash
set -e

# BenchmarkDA: Modular Master Execution Script
# Orchestrates differential abundance benchmarking with granular control

# =============================================================================
# CONFIGURATION AND USAGE
# =============================================================================

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
cd "$SCRIPT_DIR"

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
BenchmarkDA: Modular Differential Abundance Benchmarking

Usage: ./cli.sh [OPTIONS] [STEPS...]

OPTIONS:
    -h, --help              Show this help message
    -d, --datasets LIST     Comma-separated list of datasets to process
                           Available: linear,branch,cluster,covid19-pbmc,pancreas,bcr-xl,levine32
                           Default: all available datasets
    -m, --methods LIST      Comma-separated list of method types to run
                           Available: python,r,all
                           Default: all
    -e, --embeddings LIST   Comma-separated list of embeddings to use
                           Available: dm,pca,both
                           Default: both
    -s, --skip-missing      Skip missing datasets instead of prompting
    --dry-run              Show what would be executed without running

STEPS (run all if none specified):
    setup                   Setup environment and directories
    preprocess             Preprocess datasets and create embeddings
    labels                 Generate synthetic condition labels
    benchmark              Run differential abundance method benchmarks
    all                    Run all steps (default)

EXAMPLES:
    # Run complete pipeline
    ./cli.sh

    # Run only preprocessing for specific datasets
    ./cli.sh --datasets linear,branch preprocess

    # Run benchmarking for Python methods only with DM embeddings
    ./cli.sh --methods python --embeddings dm benchmark

    # Run setup and preprocessing, skip missing datasets
    ./cli.sh --skip-missing setup preprocess

    # Dry run to see what would be executed
    ./cli.sh --datasets linear --methods python --dry-run
EOF
}

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

# Default values
DATASETS=""
METHODS="all"
EMBEDDINGS="both"
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
        -e|--embeddings)
            EMBEDDINGS="$2"
            shift 2
            ;;
        -s|--skip-missing)
            SKIP_MISSING=true
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        setup|preprocess|labels|benchmark|all)
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

# Parse embeddings
if [ "$EMBEDDINGS" = "both" ]; then
    SELECTED_EMBEDDINGS=("dm" "pca")
else
    IFS=',' read -ra SELECTED_EMBEDDINGS <<< "$EMBEDDINGS"
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

    # Activate environment
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
    print_step "Preprocessing datasets with embeddings"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        local n_dm=$(get_n_dm_for_dataset "$dataset")

        print_info "Processing $dataset (DM components: $n_dm)"

        for embedding in "${SELECTED_EMBEDDINGS[@]}"; do
            if [[ "$embedding" == "dm" && $n_dm -gt 0 ]]; then
                if [ "$DRY_RUN" = true ]; then
                    echo "[DRY RUN] bash bin/dataset_preprocessing.sh \"$dataset\" X_pca \"$n_dm\" DM"
                else
                    bash bin/dataset_preprocessing.sh "$dataset" X_pca "$n_dm" DM || {
                        print_warning "DM preprocessing failed for $dataset"
                    }
                fi
            elif [[ "$embedding" == "pca" ]]; then
                if [ "$DRY_RUN" = true ]; then
                    echo "[DRY RUN] bash bin/dataset_preprocessing.sh \"$dataset\" X_pca 0 PCA"
                else
                    bash bin/dataset_preprocessing.sh "$dataset" X_pca 0 PCA || {
                        print_warning "PCA preprocessing failed for $dataset"
                    }
                fi
            fi
        done
    done

    print_success "Dataset preprocessing completed"
}

generate_labels() {
    print_step "Generating synthetic condition labels"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        local n_dm=$(get_n_dm_for_dataset "$dataset")

        print_info "Generating labels for $dataset"

        for embedding in "${SELECTED_EMBEDDINGS[@]}"; do
            if [[ "$embedding" == "dm" && $n_dm -gt 0 ]]; then
                if [ "$DRY_RUN" = true ]; then
                    echo "[DRY RUN] bash bin/modified_benchmarkda_dm_all.sh \"$dataset\" dm No DM \"$n_dm\""
                else
                    bash bin/modified_benchmarkda_dm_all.sh "$dataset" dm No DM "$n_dm" || {
                        print_warning "DM label generation failed for $dataset"
                    }
                fi
            elif [[ "$embedding" == "pca" ]]; then
                if [ "$DRY_RUN" = true ]; then
                    echo "[DRY RUN] bash bin/modified_benchmarkda.sh \"$dataset\" pca No PCA 0"
                else
                    bash bin/modified_benchmarkda.sh "$dataset" pca No PCA 0 || {
                        print_warning "PCA label generation failed for $dataset"
                    }
                fi
            fi
        done
    done

    print_success "Label generation completed"
}

run_benchmarks() {
    print_step "Running differential abundance method benchmarks"

    for dataset in "${SELECTED_DATASETS[@]}"; do
        local n_dm=$(get_n_dm_for_dataset "$dataset")

        for method_type in "${SELECTED_METHODS[@]}"; do
            for embedding in "${SELECTED_EMBEDDINGS[@]}"; do
                local analysis_layer embedding_mode n_dm_val

                if [[ "$embedding" == "dm" ]]; then
                    analysis_layer="dm"
                    embedding_mode="DM"
                    n_dm_val="$n_dm"
                else
                    analysis_layer="pca"
                    embedding_mode="PCA"
                    n_dm_val="0"
                fi

                print_info "Running $method_type methods on $dataset ($embedding_mode)"

                if [ "$DRY_RUN" = true ]; then
                    echo "[DRY RUN] python bin/run_benchmark.py --dataset \"$dataset\" --method_type \"$method_type\" --analysis_layer \"$analysis_layer\" --n_dm \"$n_dm_val\" --mode_embedding \"$embedding_mode\""
                else
                    python bin/run_benchmark.py \
                        --dataset "$dataset" \
                        --method_type "$method_type" \
                        --analysis_layer "$analysis_layer" \
                        --n_dm "$n_dm_val" \
                        --mode_embedding "$embedding_mode" \
                        --output "benchmark_scripts/${dataset}_${embedding_mode,,}_${method_type}.sh"

                    # Execute the generated script
                    local script_file="benchmark_scripts/${dataset}_${embedding_mode,,}_${method_type}.sh"
                    if [[ -f "$script_file" ]]; then
                        bash "$script_file" || {
                            print_warning "Benchmark execution failed for $dataset $method_type $embedding_mode"
                        }
                    fi
                fi
            done
        done
    done

    print_success "Benchmark execution completed"
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
    print_info "  Embeddings: ${SELECTED_EMBEDDINGS[*]}"
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