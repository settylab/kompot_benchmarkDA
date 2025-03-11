#!/bin/bash
set -e

# ====================================================================
# BenchmarkDA: Benchmarking Differential Abundance Methods in Single-Cell Data
# ====================================================================
# This master script orchestrates the entire workflow:
# 1. Setup environments for Python and R methods
# 2. Prepare datasets (download or synthetic generation)
# 3. Preprocess datasets with PCA and diffusion maps
# 4. Generate synthetic condition labels for benchmarking
# 5. Execute all DA methods on all datasets using Slurm

echo "Starting BenchmarkDA workflow..."

# Create logs directory for Slurm output
mkdir -p SlurmLog

# ====================================================================
# 1. Environment Setup
# ====================================================================
echo "Setting up computational environments..."

# Python environment (for MELD, CNA, Mellon methods)
if ! micromamba env list | grep -q "diffabundance"; then
    echo "Creating Python environment with micromamba..."
    micromamba create -n diffabundance -f differential_abundance_env_list.yml -y
fi

eval "$(micromamba shell hook --shell bash)"
micromamba activate diffabundance

# R environment (for Milo, DAseq, CyDAR methods)
echo "Setting up R environment..."
mkdir -p .Renviron
echo "R_LIBS_USER=./renv/library" > .Renviron
Rscript -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv', repos = 'https://cloud.r-project.org/'); renv::restore()"

# ====================================================================
# 2. Dataset Preparation
# ====================================================================
echo "Setting up dataset directories..."

# Create directories for all datasets
mkdir -p data/synthetic
mkdir -p data/real/bcr-xl
mkdir -p data/real/covid19-pbmc
mkdir -p data/real/levine32
mkdir -p data/real/pancreas

# Prompt for real datasets (these must be downloaded separately)
echo "IMPORTANT: Real datasets must be downloaded manually"
echo "Download from: https://drive.google.com/drive/folders/15wWFD5FMe0VdzN1pUnaUUpQ17OXkeebH"
echo "Place files in these directories:"
echo "  - data/real/bcr-xl/"
echo "  - data/real/covid19-pbmc/"
echo "  - data/real/levine32/"
echo "  - data/real/pancreas/"
echo "Press Enter to continue (or Ctrl+C to cancel)..."
read -r

# ====================================================================
# 3. Dataset Preprocessing
# ====================================================================
echo "Preprocessing datasets..."

# Process synthetic datasets
for DATASET in linear branch cluster; do
    # Create diffusion map embeddings (10 components)
    echo "Creating diffusion map for $DATASET..."
    bash bin/dataset_preprocessing.sh $DATASET X_pca 10 DM
    
    # Create PCA embeddings
    echo "Creating PCA for $DATASET..."
    bash bin/dataset_preprocessing.sh $DATASET X_pca 0 PCA
done

# Process real datasets if available
for DATASET in covid19-pbmc bcr-xl levine32 pancreas; do
    DATASET_PATH="data/real/${DATASET}/${DATASET}.h5ad"
    if [ -f "$DATASET_PATH" ]; then
        echo "Processing real dataset: $DATASET"
        
        # Real datasets use 30 diffusion components
        echo "Creating diffusion map for $DATASET..."
        bash bin/dataset_preprocessing.sh $DATASET X_pca 30 DM
        
        echo "Creating PCA for $DATASET..."
        bash bin/dataset_preprocessing.sh $DATASET X_pca 0 PCA
    else
        echo "Skipping $DATASET (file not found: $DATASET_PATH)"
    fi
done

# ====================================================================
# 4. Generate Synthetic Labels
# ====================================================================
echo "Generating condition labels for benchmarking..."

# For synthetic datasets
for DATASET in linear branch cluster; do
    # Generate labels for diffusion map embeddings
    echo "Generating labels for $DATASET with diffusion map..."
    bash bin/modified_benchmarkda_dm_all.sh $DATASET dm No DM 10
    
    # Generate labels for PCA embeddings
    echo "Generating labels for $DATASET with PCA..."
    bash bin/modified_benchmarkda.sh $DATASET pca No PCA 0
done

# For real datasets if available
for DATASET in covid19-pbmc bcr-xl levine32 pancreas; do
    if [ -f "data/real/${DATASET}/${DATASET}.h5ad" ]; then
        echo "Generating labels for $DATASET with diffusion map..."
        bash bin/modified_benchmarkda_dm_all.sh $DATASET dm No DM 30
        
        echo "Generating labels for $DATASET with PCA..."
        bash bin/modified_benchmarkda.sh $DATASET pca No PCA 0
    fi
done

# ====================================================================
# 5. Run Benchmarking with Unified Script Generator
# ====================================================================
echo "Submitting benchmark jobs to Slurm..."

# Common parameters
ITERATION=0
BALANCE="No"
TEMP_DIR="benchmark_scripts"

# Create temporary directory for generated scripts
mkdir -p $TEMP_DIR

# Function to generate and submit benchmark job
run_benchmark() {
    local dataset=$1
    local method_type=$2
    local analysis_layer=$3
    local n_dm=$4
    local mode_embedding=$5
    
    # Create descriptive script name
    local script_name="${TEMP_DIR}/benchmark_${dataset}_${mode_embedding,,}_${method_type}.sh"
    
    echo "Generating script for $dataset ($mode_embedding embedding) with $method_type methods..."
    
    # Generate the benchmark script
    python bin/run_benchmark.py \
        --dataset $dataset \
        --method_type $method_type \
        --analysis_layer $analysis_layer \
        --iteration_num $ITERATION \
        --balance $BALANCE \
        --n_dm $n_dm \
        --mode_embedding $mode_embedding \
        --output $script_name
    
    chmod +x $script_name
    
    echo "Submitting $script_name to Slurm..."
    ./$script_name
}

# Process synthetic datasets
for DATASET in linear branch cluster; do
    # Python methods
    run_benchmark $DATASET "python" "dm" 10 "DM"
    run_benchmark $DATASET "python" "pca" 0 "PCA"
    
    # R methods
    run_benchmark $DATASET "r" "dm" 10 "DM"
    run_benchmark $DATASET "r" "pca" 0 "PCA"
done

# Process real datasets if available
for DATASET in covid19-pbmc bcr-xl levine32 pancreas; do
    if [ -f "data/real/${DATASET}/${DATASET}.h5ad" ]; then
        # Python methods
        run_benchmark $DATASET "python" "dm" 30 "DM"
        run_benchmark $DATASET "python" "pca" 0 "PCA"
        
        # R methods
        run_benchmark $DATASET "r" "dm" 30 "DM"
        run_benchmark $DATASET "r" "pca" 0 "PCA"
    fi
done

echo "====================================================================================="
echo "Benchmark jobs submitted to Slurm"
echo "====================================================================================="
echo "EXPLANATION:"
echo "1. Runtime-generated scripts are stored in the '$TEMP_DIR' directory"
echo "2. Each script submits an array job to Slurm with all parameter combinations"
echo "3. Job logs will be saved to 'SlurmLog/[jobname]_[node]_[jobid]_[taskid].out'"
echo "4. Results will be saved to 'benchmark/[embedding]/[dataset]/[parameters]/'"
echo ""
echo "Check job status with: squeue -u $USER"
echo "====================================================================================="
