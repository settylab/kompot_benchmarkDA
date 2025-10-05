#!/bin/bash

# Note: Environment activation is handled by the calling script (cli.sh)
# This script relies on MAMBA_EXE being set and uses `mamba run` for all Python commands

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"
if [ -z "${root+x}" ]; then
    export root="$(readlink -f "$script_dir/..")"
fi
# Stay in project root (no cd to python_method - that directory doesn't exist) 

#data_dir="$root/data"
#out_dir="$root/benchmark_python"

#echo "data_dir is $data_dir" 

## Run real data ##

# Parse arguments
data_id=$1
embedding_layer=$2
n_dm=$3
mode_embedding=$4

# Determine dataset type (synthetic vs real)
if [[ "$data_id" == "linear" || "$data_id" == "branch" || "$data_id" == "cluster" ]]; then
    data_dir=${root}/data/synthetic/$data_id
    data_type="synthetic"
else
    data_dir=${root}/data/real/$data_id
    data_type="real"
fi

# Print step header
echo ""
echo -e "${BLUE}===================================================${NC}"
if [ "$mode_embedding" = "PCA" ]; then
    echo -e "${BLUE}Step 1/2: Creating X_pca embeddings${NC}"
    echo -e "${BLUE}Dataset: ${data_id} (${data_type})${NC}"
else
    echo -e "${BLUE}Step 2/2: Computing DM from X_pca (${n_dm} components)${NC}"
    echo -e "${BLUE}Dataset: ${data_id} (${data_type})${NC}"
fi
echo -e "${BLUE}===================================================${NC}"

# Get pop_col from config using Python
pop_col=$(python -c "
import sys
sys.path.insert(0, '${root}')
from config.dataset_config import DATASET_CONFIGS
config = DATASET_CONFIGS.get('$data_id', {})
print(config.get('pop_col', 'celltype'))
" 2>/dev/null)

# Fallback if Python fails
if [ -z "$pop_col" ]; then
    echo -e "${YELLOW}[WARN] Could not read pop_col from config, using default 'celltype'${NC}"
    pop_col="celltype"
fi

# Set file paths
data_file=${data_dir}/${data_id}.h5ad
out_dir=${data_dir}/${data_id}.h5ad.tmp  # Write to temp file first

echo -e "${BLUE}[INFO] Input file: ${data_file}${NC}"
echo -e "${BLUE}[INFO] Working directory: ${data_dir}${NC}"
echo ""

# Run preprocessing pipeline
echo -e "${BLUE}[STEP] Running preprocessing pipeline...${NC}"
if [ -n "${MAMBA_EXE}" ]; then
    ${MAMBA_EXE} run -n benchmarkda python ${root}/lib/data_preprocessing_pipeline.py \
        --file_path "${data_file}" \
        --embedding_layer "${embedding_layer}" \
        --n_dm "${n_dm}" \
        --pop_col "${pop_col}" \
        --mode_embedding "${mode_embedding}" \
        --output_dir "${out_dir}" 2>&1 | grep -v "^During startup" | grep -v "^package 'colorout'"
else
    python ${root}/lib/data_preprocessing_pipeline.py \
        --file_path "${data_file}" \
        --embedding_layer "${embedding_layer}" \
        --n_dm "${n_dm}" \
        --pop_col "${pop_col}" \
        --mode_embedding "${mode_embedding}" \
        --output_dir "${out_dir}" 2>&1 | grep -v "^During startup" | grep -v "^package 'colorout'"
fi

echo ""

# Skip RDS conversion - R methods read H5AD directly via R anndata package
# (scripts/run_DA.r line 75: adata <- read_h5ad(args$file_path))
echo -e "${BLUE}[INFO] Skipping RDS conversion (R methods use H5AD directly)${NC}"

# Move temp file to main file (overwrite)
if [ -f "${out_dir}" ]; then
    echo -e "${BLUE}[STEP] Updating main file: ${data_file}${NC}"
    mv "${out_dir}" "${data_file}"
fi

echo ""
if [ "$mode_embedding" = "DM" ]; then
    echo -e "${GREEN}[DONE] Preprocessing complete: X_pca ($(cat ${data_file} 2>/dev/null | wc -c | awk '{print int($1/1024/1024)}')MB) + DM_EigenVectors (${n_dm} components)${NC}"
else
    echo -e "${GREEN}[DONE] X_pca embeddings created${NC}"
fi
echo ""