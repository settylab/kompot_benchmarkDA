#!/bin/bash
set -e

# Comprehensive script for benchmarkDA setup and execution
# This script automates the entire benchmark process from environment setup to result generation
# Based on examining the project structure and individual scripts in bin/

echo "Starting benchmarkDA setup and execution..."

# --------------------------------
# 1. Environment Setup
# --------------------------------
# This section creates and activates the necessary Python and R environments
# Python methods (MELD, CNA, Mellon) need the micromamba environment
# R methods (Milo, DAseq, CyDAR) need the renv environment
echo "Setting up environment..."

# Create micromamba environment for Python dependencies
# The differential_abundance_env_list.yml file contains all required Python packages
micromamba create -n diffabundance -f differential_abundance_env_list.yml -y
eval "$(micromamba shell hook --shell bash)"
micromamba activate diffabundance

# Set up R environment using renv
# renv.lock contains all required R packages including Milo, DAseq, and CyDAR
mkdir -p .Renviron
echo "R_LIBS_USER=./renv/library" > .Renviron
Rscript -e "if (!requireNamespace('renv', quietly = TRUE)) install.packages('renv', repos = 'https://cloud.r-project.org/'); renv::restore()"

# --------------------------------
# 2. Download datasets
# --------------------------------
# This section creates the data directory and prompts for downloading real datasets
# Real datasets (bcr-xl, covid19-pbmc, levine32, and pancreas) need to be downloaded manually
# From examining R_method_evaluation_dm.sh and modified_benchmarkda.sh:
# - bcr-xl should be in data/real/bcr-xl/ directory
# - covid19-pbmc should be in data/real/covid19-pbmc/ directory
# - levine32 should be in data/real/levine32/ directory
# - pancreas should be in data/real/pancreas/ directory
echo "Setting up data directories..."

# Create data directory if it doesn't exist
mkdir -p data
mkdir -p data/real/bcr-xl
mkdir -p data/real/covid19-pbmc
mkdir -p data/real/levine32
mkdir -p data/real/pancreas

# For real datasets, you need to download them manually from Google Drive
echo "Please download real datasets from: https://drive.google.com/drive/folders/15wWFD5FMe0VdzN1pUnaUUpQ17OXkeebH"
echo "and place them in the correct directories:"
echo "- bcr-xl dataset in data/real/bcr-xl/"
echo "- covid19-pbmc dataset in data/real/covid19-pbmc/"
echo "- levine32 dataset in data/real/levine32/"
echo "- pancreas dataset in data/real/pancreas/"
echo "Press Enter when done or if you want to skip real datasets and only use synthetic ones"
read -r

# --------------------------------
# 3. Dataset preprocessing
# --------------------------------
# This section converts raw datasets into preprocessed formats with both
# diffusion map (DM) and principal component analysis (PCA) representations
# From examining dataset_preprocessing.sh and the benchmark scripts:
# - Synthetic datasets use 10 diffusion components
# - Real datasets use 30 diffusion components
# - PCA uses all components (n_dm = 0)
# - X_pca is the dataset layer containing PCA embeddings
echo "Starting dataset preprocessing..."

# Process synthetic datasets with diffusion map and PCA
# Looking at bin/modified_benchmarkda.sh, synthetic datasets are in data/synthetic/
for DATASET in linear branch cluster
do
    echo "Preprocessing $DATASET dataset with diffusion map..."
    # Arguments: dataset_name, embedding_layer, n_components, embedding_type
    bash bin/dataset_preprocessing.sh $DATASET X_pca 10 DM
    
    echo "Preprocessing $DATASET dataset with PCA..."
    # For PCA, n_components=0 means use all components
    bash bin/dataset_preprocessing.sh $DATASET X_pca 0 PCA
done

# Process real datasets if they exist
# Looking at bin/R_method_evaluation_dm.sh, real datasets require specific paths
for DATASET in covid19-pbmc bcr-xl levine32 pancreas
do
    # Check correct path based on script examination
    if [ -f "data/real/${DATASET}/${DATASET}.h5ad" ]; then
        echo "Preprocessing $DATASET dataset with diffusion map..."
        # Real datasets use 30 diffusion components
        bash bin/dataset_preprocessing.sh $DATASET X_pca 30 DM
        
        echo "Preprocessing $DATASET dataset with PCA..."
        bash bin/dataset_preprocessing.sh $DATASET X_pca 0 PCA
    else
        echo "Skipping $DATASET - dataset file not found at data/real/${DATASET}/${DATASET}.h5ad"
    fi
done

# --------------------------------
# 4. Generate synthetic labels
# --------------------------------
# This section creates synthetic differential abundance labels for benchmarking
# From examining modified_benchmarkda.sh and modified_benchmarkda_dm_all.sh:
# - These scripts generate synthetic condition labels using different parameters
# - "No" argument refers to BALANCE parameter (no balancing between conditions)
# - The scripts submit Slurm jobs for each parameter combination
# - Specifically, they run generate_bm_data.py to create synthetic labels
echo "Generating synthetic labels..."

# For synthetic datasets with diffusion map
for DATASET in linear branch cluster
do
    echo "Generating synthetic labels for $DATASET with diffusion map..."
    # Arguments: dataset_name, analysis_layer, balance_bool, mode_embedding, n_dm
    bash bin/modified_benchmarkda_dm_all.sh $DATASET dm No DM 10
    
    echo "Generating synthetic labels for $DATASET with PCA..."
    bash bin/modified_benchmarkda.sh $DATASET pca No PCA 0
done

# For real datasets with diffusion map if they exist
for DATASET in covid19-pbmc bcr-xl levine32 pancreas
do
    if [ -f "data/real/${DATASET}/${DATASET}.h5ad" ]; then
        echo "Generating synthetic labels for $DATASET with diffusion map..."
        # Real datasets use 30 diffusion components
        bash bin/modified_benchmarkda_dm_all.sh $DATASET dm No DM 30
        
        echo "Generating synthetic labels for $DATASET with PCA..."
        bash bin/modified_benchmarkda.sh $DATASET pca No PCA 0
    else
        echo "Skipping $DATASET - dataset file not found"
    fi
done

# --------------------------------
# 5. Run benchmarking for Python methods
# --------------------------------
# This section runs the Python-based differential abundance methods
# From examining bm_syn_python_*.sh scripts, these run:
# - MELD (meld_bm.py)
# - CNA (CNA_bm.py)
# - Mellon (Mellon_bm.py)
echo "Running benchmarking for Python methods..."

# Parameters for all benchmarking runs
# These parameters were determined by examining the benchmark scripts:
# - ITERATION=0: Only run one iteration per parameter set
# - BALANCE="No": Do not force balanced conditions
# - MELLON_METHOD="fractal": Use fractal-based algorithm for Mellon method
# - MELLON_NORM="No": Do not normalize Mellon results
# - MELLON_CORRECTION="No": Do not use correction in Mellon
# - MELLON_SYNC="Yes": Use synchronization in Mellon
ITERATION=0
BALANCE="No"
MELLON_METHOD="fractal" 
MELLON_NORM="No"
MELLON_CORRECTION="No"
MELLON_SYNC="Yes"

# For synthetic datasets with diffusion map
for DATASET in linear branch cluster
do
    echo "Benchmarking Python methods on $DATASET with diffusion map..."
    # Each benchmark script submits Slurm jobs with these parameters
    bash bin/bm_syn_python_synthetic_dm.sh $DATASET dm $ITERATION $BALANCE $MELLON_METHOD $MELLON_NORM $MELLON_CORRECTION $MELLON_SYNC 10 DM
    
    echo "Benchmarking Python methods on $DATASET with PCA..."
    bash bin/bm_syn_python_synthetic_pca_final.sh $DATASET pca $ITERATION $BALANCE $MELLON_METHOD $MELLON_NORM $MELLON_CORRECTION $MELLON_SYNC 0 PCA
done

# For real datasets if they exist
for DATASET in covid19-pbmc bcr-xl levine32 pancreas
do
    if [ -f "data/real/${DATASET}/${DATASET}.h5ad" ]; then
        echo "Benchmarking Python methods on $DATASET with diffusion map..."
        bash bin/bm_syn_python_real_dm.sh $DATASET dm $ITERATION $BALANCE $MELLON_METHOD $MELLON_NORM $MELLON_CORRECTION $MELLON_SYNC 30 DM
        
        echo "Benchmarking Python methods on $DATASET with PCA..."
        bash bin/bm_syn_python_real_pca_final.sh $DATASET pca $ITERATION $BALANCE $MELLON_METHOD $MELLON_NORM $MELLON_CORRECTION $MELLON_SYNC 0 PCA
    else
        echo "Skipping $DATASET - dataset file not found"
    fi
done

# --------------------------------
# 6. Run benchmarking for R methods
# --------------------------------
# This section runs the R-based differential abundance methods
# From examining R_method_evaluation_*.sh and run_DA.r:
# - Milo: Neighborhood-based DA method
# - DAseq: Differential abundance of cell states
# - CyDAR: Cytometry-based DA method
# - Louvain: Clustering-based DA analysis
echo "Running benchmarking for R methods..."

# For synthetic datasets with diffusion map
for DATASET in linear branch cluster
do
    echo "Benchmarking R methods on $DATASET with diffusion map..."
    # Use same parameters as Python methods for consistency
    # These scripts call run_DA.r with each method (milo, daseq, cydar, louvain)
    # Each evaluation script submits Slurm jobs
    bash bin/R_method_evaluation_dm.sh $DATASET dm $ITERATION $BALANCE $MELLON_METHOD $MELLON_NORM $MELLON_CORRECTION $MELLON_SYNC 10 DM
    
    echo "Benchmarking R methods on $DATASET with PCA..."
    bash bin/R_method_evaluation_pca.sh $DATASET pca $ITERATION $BALANCE $MELLON_METHOD $MELLON_NORM $MELLON_CORRECTION $MELLON_SYNC 0 PCA
done

# For real datasets if they exist
for DATASET in covid19-pbmc bcr-xl levine32 pancreas
do
    if [ -f "data/real/${DATASET}/${DATASET}.h5ad" ]; then
        echo "Benchmarking R methods on $DATASET with diffusion map..."
        bash bin/R_method_evaluation_dm.sh $DATASET dm $ITERATION $BALANCE $MELLON_METHOD $MELLON_NORM $MELLON_CORRECTION $MELLON_SYNC 30 DM
        
        echo "Benchmarking R methods on $DATASET with PCA..."
        bash bin/R_method_evaluation_pca.sh $DATASET pca $ITERATION $BALANCE $MELLON_METHOD $MELLON_NORM $MELLON_CORRECTION $MELLON_SYNC 0 PCA
    else
        echo "Skipping $DATASET - dataset file not found"
    fi
done

# All benchmarking scripts submit Slurm jobs with the sbatch command
# Jobs are logged to the SlurmLog/ directory
# Results are saved to benchmark/ directory for each dataset and method

echo "Benchmark execution completed!"
echo "Results can be found in the benchmark/ directory"
