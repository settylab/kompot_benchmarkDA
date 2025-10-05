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

# Configuration - derived from dataset_config.py
# Synthetic datasets are those listed in get_data_file_path as synthetic
SYNTHETIC_DATASETS=(linear branch cluster)
# Get all datasets from config and filter out synthetic ones for real datasets
ALL_DATASETS=($(sed -n '/^DATASET_CONFIGS = {/,/^}/p' config/dataset_config.py | grep '^\s*"[^"]*": {$' | sed 's/.*"\([^"]*\)".*/\1/'))
REAL_DATASETS=()
for ds in "${ALL_DATASETS[@]}"; do
    if ! echo "${SYNTHETIC_DATASETS[@]}" | grep -qw "$ds"; then
        REAL_DATASETS+=("$ds")
    fi
done

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
    # Get available methods and datasets dynamically from config files (no Python needed)
    local python_methods=$(sed -n '/PYTHON_METHODS = {/,/^R_METHODS/p' config/method_config.py | grep '^\s*"[^"]*": {$' | grep -v '"params"' | sed 's/.*"\([^"]*\)".*/\1/' | tr '\n' ',' | sed 's/,$//')
    local r_methods=$(sed -n '/^R_METHODS = {/,/^}/p' config/method_config.py | grep '^\s*"[^"]*": {$' | grep -v '"params"' | sed 's/.*"\([^"]*\)".*/\1/' | tr '\n' ',' | sed 's/,$//')
    local all_datasets=$(echo "${ALL_DATASETS[@]}" | tr ' ' ',')

    cat << EOF
BenchmarkDA: Differential Abundance Benchmarking Pipeline

USAGE: ./cli.sh [OPTIONS] [STEPS...]

OPTIONS:
    -h, --help              Show this help message
    -d, --datasets LIST     Datasets: $all_datasets (from config/dataset_config.py)
    -m, --methods LIST      Methods: python,r,all OR specific method names (default: all)
    -s, --skip-missing      Skip missing datasets without prompting
    --only-missing          Generate only missing labels (skip existing combinations)
    --slurm                Submit benchmarks as SLURM array jobs
    --sbatch-options OPTS  Forward SLURM options to sbatch (e.g., "--partition=gpu --gres=gpu:1")
    --dry-run              Show commands without executing

FILTERING OPTIONS (for labels step):
    --populations LIST      Filter by populations (e.g., M2,M8)
    --seeds LIST           Filter by seeds (e.g., 43,44)
    --enrichments LIST     Filter by enrichment values (e.g., 0.75,0.95)
    --batch-sds LIST       Filter by batch standard deviations (e.g., 0.75,1.25)

STEPS:
    setup                   Setup environment and directories
    convert-data           Convert RDS files to H5AD format
    preprocess             Create PCA embeddings and DM from PCA+batch
    labels                 Generate synthetic labels
    benchmark              Run DA method benchmarks
    test                   Run test suite
    status                 Show completion status
    all                    Run all steps (default)

STATUS COMMANDS:
    ./cli.sh status                              All datasets
    ./cli.sh --datasets linear,branch status     Specific datasets

STATUS OUTPUT:
    [Complete] = Finished   [Partial] = In progress   [Missing] = Not started

PYTHON METHODS:
    ${python_methods}

R METHODS:
    ${r_methods}

EXAMPLES:
    ./cli.sh                                         Complete pipeline (local)
    ./cli.sh convert-data                            Convert all RDS files to H5AD
    ./cli.sh test                                    Run test suite
    ./cli.sh labels --slurm                          Submit label generation as SLURM array (splits by population)
    ./cli.sh --datasets bcr-xl labels --slurm        Submit labels for specific dataset(s)
    ./cli.sh benchmark --slurm                       Submit ALL benchmarks to SLURM
    ./cli.sh --methods milo benchmark --slurm --sbatch-options "--partition=largenode --mem=64G"    Submit with custom SLURM options
    ./cli.sh --methods milo benchmark --slurm        Submit only Milo method
    ./cli.sh --datasets branch --methods milo benchmark --slurm    Small test: branch + milo
    ./cli.sh --datasets branch --populations M2,M8 --seeds 43,44 --enrichments 0.75,0.95 --batch-sds 0.75,1.25,1.5 labels    Generate labels for specific combinations
    ./cli.sh status                                  Check progress
    ./cli.sh --datasets linear preprocess            Preprocess one dataset
    ./cli.sh --methods python benchmark              Python methods only
    ./cli.sh --methods r benchmark --slurm           Submit only R methods to SLURM
    ./cli.sh --dry-run labels --slurm                Show what SLURM scripts would be generated
EOF
}

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

# Default values
DATASETS=""
METHODS="all"
SKIP_MISSING=false
ONLY_MISSING=false
DRY_RUN=false
STEPS=()
FILTER_POPULATIONS=""
FILTER_SEEDS=""
FILTER_ENRICHMENTS=""
FILTER_BATCH_SDS=""
SBATCH_OPTIONS=""

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
        --only-missing)
            ONLY_MISSING=true
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
        --populations)
            FILTER_POPULATIONS="$2"
            shift 2
            ;;
        --seeds)
            FILTER_SEEDS="$2"
            shift 2
            ;;
        --enrichments)
            FILTER_ENRICHMENTS="$2"
            shift 2
            ;;
        --batch-sds)
            FILTER_BATCH_SDS="$2"
            shift 2
            ;;
        --sbatch-options)
            SBATCH_OPTIONS="$2"
            shift 2
            ;;
        setup|convert-data|preprocess|labels|benchmark|test|status|all)
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

# Get available methods from config file directly (no Python needed)
AVAILABLE_PYTHON_METHODS=($(sed -n '/PYTHON_METHODS = {/,/^R_METHODS/p' config/method_config.py | grep '^\s*"[^"]*": {$' | grep -v '"params"' | sed 's/.*"\([^"]*\)".*/\1/'))
AVAILABLE_R_METHODS=($(sed -n '/^R_METHODS = {/,/^}/p' config/method_config.py | grep '^\s*"[^"]*": {$' | grep -v '"params"' | sed 's/.*"\([^"]*\)".*/\1/'))

# Parse methods - support both method types (python/r) and individual method names (milo, kompot, etc)
SPECIFIC_METHODS=()
if [ "$METHODS" = "all" ]; then
    SELECTED_METHODS=("python" "r")
else
    IFS=',' read -ra METHOD_LIST <<< "$METHODS"
    SELECTED_METHODS=()
    for method in "${METHOD_LIST[@]}"; do
        case "$method" in
            python|r)
                SELECTED_METHODS+=("$method")
                ;;
            *)
                # Check if it's a valid Python method
                if [[ " ${AVAILABLE_PYTHON_METHODS[@]} " =~ " $method " ]]; then
                    if [[ ! " ${SELECTED_METHODS[@]} " =~ " python " ]]; then
                        SELECTED_METHODS+=("python")
                    fi
                    SPECIFIC_METHODS+=("$method")
                # Check if it's a valid R method
                elif [[ " ${AVAILABLE_R_METHODS[@]} " =~ " $method " ]]; then
                    if [[ ! " ${SELECTED_METHODS[@]} " =~ " r " ]]; then
                        SELECTED_METHODS+=("r")
                    fi
                    SPECIFIC_METHODS+=("$method")
                else
                    print_error "Unknown method: $method"
                    print_error "Available Python methods: ${AVAILABLE_PYTHON_METHODS[*]}"
                    print_error "Available R methods: ${AVAILABLE_R_METHODS[*]}"
                    exit 1
                fi
                ;;
        esac
    done
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
    # Create directories for all datasets from config
    mkdir -p data/synthetic data/real SlurmLog benchmark_scripts
    for ds in "${REAL_DATASETS[@]}"; do
        mkdir -p "data/real/$ds"
    done

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

convert_data() {
    print_step "Converting RDS files to H5AD format"

    # Conversion needs environment active
    if [ "$DRY_RUN" != true ]; then
        # Load environment utilities
        if [[ -f "bin/environment_utils.sh" ]]; then
            source "bin/environment_utils.sh"
        fi

        # Activate environment if not already active
        if [ -z "$CONDA_PREFIX" ]; then
            print_info "Activating environment for data conversion"
            if ! activate_benchmarkda_environment; then
                print_error "Failed to activate environment"
                return 1
            fi
        fi
    fi

    local converted_count=0
    local skipped_count=0

    for dataset in "${SELECTED_DATASETS[@]}"; do
        # Skip synthetic datasets
        if [[ " ${SYNTHETIC_DATASETS[*]} " =~ " ${dataset} " ]]; then
            continue
        fi

        local rds_file="data/real/$dataset/$dataset.rds"
        local h5ad_file="data/real/$dataset/$dataset.h5ad"

        if [[ ! -f "$rds_file" ]]; then
            print_warning "RDS file not found: $rds_file"
            continue
        fi

        if [[ -f "$h5ad_file" ]]; then
            print_info "H5AD already exists for $dataset, skipping"
            ((skipped_count++))
            continue
        fi

        print_info "Converting $dataset: RDS → H5AD"

        if [ "$DRY_RUN" = true ]; then
            echo "[DRY RUN] python bin/convert_rds_to_h5ad.py --input \"$rds_file\" --output \"$h5ad_file\""
        else
            python bin/convert_rds_to_h5ad.py --input "$rds_file" --output "$h5ad_file" || {
                print_warning "Conversion failed for $dataset"
                continue
            }
            ((converted_count++))
        fi
    done

    if [ "$DRY_RUN" != true ]; then
        print_success "Data conversion completed: $converted_count converted, $skipped_count skipped"
    fi
}

preprocess_datasets() {
    print_step "Preprocessing datasets"

    # Preprocessing needs environment active
    if [ "$DRY_RUN" != true ]; then
        # Load environment utilities
        if [[ -f "bin/environment_utils.sh" ]]; then
            source "bin/environment_utils.sh"
        fi

        # Activate environment if not already active
        if [ -z "$CONDA_PREFIX" ]; then
            print_info "Activating environment for preprocessing"
            if ! activate_benchmarkda_environment; then
                print_error "Failed to activate environment"
                return 1
            fi
        fi
    fi

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

    # Check if SLURM submission is requested
    if [ "$USE_SLURM" = true ]; then
        generate_labels_slurm
        return
    fi

    # Label generation needs environment active
    if [ "$DRY_RUN" != true ]; then
        # Load environment utilities
        if [[ -f "bin/environment_utils.sh" ]]; then
            source "bin/environment_utils.sh"
        fi

        # Activate environment if not already active
        if [ -z "$CONDA_PREFIX" ]; then
            print_info "Activating environment for label generation"
            if ! activate_benchmarkda_environment; then
                print_error "Failed to activate environment"
                return 1
            fi
        fi
    fi

    # Build filter arguments (all optional)
    local filter_args=""
    [ -n "$FILTER_POPULATIONS" ] && filter_args="$filter_args --populations $FILTER_POPULATIONS"
    [ -n "$FILTER_SEEDS" ] && filter_args="$filter_args --seeds $FILTER_SEEDS"
    [ -n "$FILTER_ENRICHMENTS" ] && filter_args="$filter_args --enrichments $FILTER_ENRICHMENTS"
    [ -n "$FILTER_BATCH_SDS" ] && filter_args="$filter_args --batch-sds $FILTER_BATCH_SDS"

    if [ -n "$filter_args" ]; then
        print_info "Using filtered label generation:"
        [ -n "$FILTER_POPULATIONS" ] && print_info "  Populations: $FILTER_POPULATIONS"
        [ -n "$FILTER_SEEDS" ] && print_info "  Seeds: $FILTER_SEEDS"
        [ -n "$FILTER_ENRICHMENTS" ] && print_info "  Enrichments: $FILTER_ENRICHMENTS"
        [ -n "$FILTER_BATCH_SDS" ] && print_info "  Batch SDs: $FILTER_BATCH_SDS"
    fi

    for dataset in "${SELECTED_DATASETS[@]}"; do
        print_info "Generating labels for $dataset"

        # Add skip-existing flag if only-missing is set
        # (CLI uses --only-missing, but Python script uses --skip-existing)
        local skip_flag=""
        if [ "$ONLY_MISSING" = true ]; then
            skip_flag="--skip-existing"
        fi

        if [ "$DRY_RUN" = true ]; then
            echo "[DRY RUN] python bin/generate_labels.py --dataset \"$dataset\" $filter_args $skip_flag"
        else
            # Use unified utility (environment already activated above)
            python bin/generate_labels.py --dataset "$dataset" $filter_args $skip_flag || {
                print_warning "Label generation failed for $dataset"
            }
        fi
    done

    print_success "Label generation completed"
}

generate_labels_slurm() {
    print_step "Generating SLURM array job scripts for label generation"

    mkdir -p benchmark_scripts SlurmLog

    # Activate environment for getting populations from config
    if [ "$DRY_RUN" != true ]; then
        if [[ -f "bin/environment_utils.sh" ]]; then
            source "bin/environment_utils.sh"
        fi
        if [ -z "$CONDA_PREFIX" ]; then
            if ! activate_benchmarkda_environment; then
                print_error "Failed to activate environment"
                return 1
            fi
        fi
    fi

    for dataset in "${SELECTED_DATASETS[@]}"; do
        # Get populations for this dataset
        local pops_array=()
        if [ -n "$FILTER_POPULATIONS" ]; then
            # Use filtered populations
            IFS=',' read -ra pops_array <<< "$FILTER_POPULATIONS"
        else
            # Get all populations from config
            pops_array=($(python -c "
import sys
sys.path.append('.')
from config.dataset_config import DATASET_CONFIGS
pops = DATASET_CONFIGS.get('$dataset', {}).get('pops', [])
print(' '.join(pops))
" 2>/dev/null))
        fi

        if [ ${#pops_array[@]} -eq 0 ]; then
            print_warning "No populations found for $dataset, skipping"
            continue
        fi

        # Build filter arguments (excluding populations since we're splitting on that)
        local filter_args=""
        [ -n "$FILTER_SEEDS" ] && filter_args="$filter_args --seeds $FILTER_SEEDS"
        [ -n "$FILTER_ENRICHMENTS" ] && filter_args="$filter_args --enrichments $FILTER_ENRICHMENTS"
        [ -n "$FILTER_BATCH_SDS" ] && filter_args="$filter_args --batch-sds $FILTER_BATCH_SDS"

        # Add only-missing flag if set
        local skip_flag=""
        if [ "$ONLY_MISSING" = true ]; then
            skip_flag="--only-missing"
        fi

        local script_path="benchmark_scripts/slurm_labels_${dataset}.sh"
        local num_pops=${#pops_array[@]}
        local max_array_idx=$((num_pops - 1))

        print_info "Creating SLURM array job for $dataset ($num_pops populations)"

        # Create the SLURM script
        cat > "$script_path" << EOF
#!/bin/bash
#SBATCH --job-name=labels_${dataset}
#SBATCH --array=0-${max_array_idx}
#SBATCH --cpus-per-task=8
#SBATCH --time=2-00:00:00
#SBATCH --output=SlurmLog/%x_%A_%a.out
#SBATCH --error=SlurmLog/%x_%A_%a.err

# BenchmarkDA Label Generation - SLURM Array Job
# Generated: $(date)
# Dataset: ${dataset}
# Populations: ${pops_array[*]}

# Define populations array
pops=(${pops_array[@]@Q})

# Get population for this task
pop=\${pops[\$SLURM_ARRAY_TASK_ID]}

echo "=========================================="
echo "SLURM Array Job: Label Generation"
echo "=========================================="
echo "Job ID: \$SLURM_JOB_ID"
echo "Array Task ID: \$SLURM_ARRAY_TASK_ID"
echo "Dataset: ${dataset}"
echo "Population: \$pop"
echo "=========================================="
echo ""

# Change to project directory
cd "${SCRIPT_DIR}"

# Run label generation for this population
./cli.sh --datasets ${dataset} --populations "\$pop" $filter_args $skip_flag labels

echo ""
echo "=========================================="
echo "Task completed: \$pop"
echo "=========================================="
EOF

        chmod +x "$script_path"

        if [ "$DRY_RUN" = true ]; then
            echo "[DRY RUN] Would create: $script_path"
            echo "[DRY RUN] Would submit: sbatch $SBATCH_OPTIONS $script_path"
        else
            print_success "Created SLURM script: $script_path"

            # Submit the job
            local sbatch_cmd="sbatch --export=ALL"
            if [ -n "$SBATCH_OPTIONS" ]; then
                sbatch_cmd="$sbatch_cmd $SBATCH_OPTIONS"
                print_info "Using custom SLURM options: $SBATCH_OPTIONS"
            fi

            local job_id=$($sbatch_cmd "$script_path" | awk '{print $NF}')
            print_success "Submitted SLURM array job $job_id for $dataset (${num_pops} populations)"
        fi
    done

    print_success "SLURM array job submission completed"
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

            # Build methods argument if specific methods were requested (do this before dry-run)
            local methods_arg=""
            if [ ${#SPECIFIC_METHODS[@]} -gt 0 ]; then
                # Filter specific methods for this method_type using config-derived arrays
                local type_methods=()
                for method in "${SPECIFIC_METHODS[@]}"; do
                    # Check if method is in R methods array
                    if echo "${AVAILABLE_R_METHODS[@]}" | grep -qw "$method"; then
                        if [ "$method_type" = "r" ]; then
                            type_methods+=("$method")
                        fi
                    # Otherwise assume it's a Python method
                    elif echo "${AVAILABLE_PYTHON_METHODS[@]}" | grep -qw "$method"; then
                        if [ "$method_type" = "python" ]; then
                            type_methods+=("$method")
                        fi
                    fi
                done

                if [ ${#type_methods[@]} -gt 0 ]; then
                    methods_arg="--methods ${type_methods[@]}"
                fi
            fi

            # Add skip-existing flag if only-missing is set (do this before dry-run)
            local skip_flag=""
            if [ "$ONLY_MISSING" = true ]; then
                skip_flag="--skip-existing"
            fi

            # Build filter arguments (all optional)
            local filter_args=""
            [ -n "$FILTER_POPULATIONS" ] && filter_args="$filter_args --populations $FILTER_POPULATIONS"
            [ -n "$FILTER_SEEDS" ] && filter_args="$filter_args --seeds $FILTER_SEEDS"
            [ -n "$FILTER_ENRICHMENTS" ] && filter_args="$filter_args --enrichments $FILTER_ENRICHMENTS"
            [ -n "$FILTER_BATCH_SDS" ] && filter_args="$filter_args --batch-sds $FILTER_BATCH_SDS"

            if [ "$DRY_RUN" = true ]; then
                echo "[DRY RUN] python bin/direct_benchmark.py --dataset \"$dataset\" --method_type \"$method_type\" --n_dm \"$n_dm\" $methods_arg $skip_flag $slurm_flag $filter_args"
            else
                if [ "$USE_SLURM" = true ]; then
                    print_info "Submitting $method_type methods for $dataset to SLURM"
                else
                    print_info "Executing benchmarks"
                fi

                # Activate environment for both SLURM and local runs
                # Load GCC module if needed for R methods (before mamba activation)
                if [ "$method_type" = "r" ]; then
                    load_gcc_module_if_needed
                fi

                # Activate environment
                if ! activate_benchmarkda_environment; then
                    print_error "Failed to activate environment"
                    continue
                fi

                # Generate SLURM script (uses activated environment)
                python bin/direct_benchmark.py \
                    --dataset "$dataset" \
                    --method_type "$method_type" \
                    --n_dm "$n_dm" \
                    $methods_arg \
                    $skip_flag \
                    $slurm_flag \
                    $filter_args || {
                    print_warning "Direct benchmark execution failed for $dataset $method_type"
                }

                # If using SLURM, submit the generated script
                if [ "$USE_SLURM" = true ] && [ "$DRY_RUN" != true ]; then
                    local script_path="benchmark_scripts/slurm_benchmark_${dataset}_${method_type}.sh"
                    if [ -f "$script_path" ]; then
                        # Build sbatch command with optional custom options
                        local sbatch_cmd="sbatch --export=ALL"
                        if [ -n "$SBATCH_OPTIONS" ]; then
                            sbatch_cmd="$sbatch_cmd $SBATCH_OPTIONS"
                            print_info "Using custom SLURM options: $SBATCH_OPTIONS"
                        fi
                        local job_id=$($sbatch_cmd "$script_path" | awk '{print $NF}')
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

run_tests() {
    print_step "Running test suite"

    # Testing needs environment active
    if [ "$DRY_RUN" != true ]; then
        # Load environment utilities
        if [[ -f "bin/environment_utils.sh" ]]; then
            source "bin/environment_utils.sh"
        fi

        # Activate environment if not already active
        if [ -z "$CONDA_PREFIX" ]; then
            print_info "Activating environment for testing"
            if ! activate_benchmarkda_environment; then
                print_error "Failed to activate environment"
                return 1
            fi
        fi
    fi

    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] python bin/run_tests.py"
    else
        print_info "Running all tests"
        python bin/run_tests.py || {
            print_error "Tests failed"
            return 1
        }
    fi

    print_success "All tests passed"
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
            convert-data)
                convert_data
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
            test)
                run_tests
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
    echo "# Run specific method on specific dataset:"
    echo "./cli.sh --datasets linear --methods kompot benchmark"
    echo ""
    echo "# Run specific components:"
    echo "./cli.sh --datasets linear --methods python preprocess benchmark"

    print_header "All Done!"
}

main "$@"