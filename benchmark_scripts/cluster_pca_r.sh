#!/bin/bash

# Note: This script is designed to be user-agnostic and does not use the module system

# Set slurm parameters
time=1-00:00:00
partition=campus-new

script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"
if [ -z "${root+x}" ]; then
    export root="$(readlink -f "$script_dir/..")"
fi

# Load environment detection utilities
source "${root}/bin/environment_utils.sh"

cd ${root}/scripts
echo "Working directory: $root/scripts"

# Parameters from command line
data_id="cluster"
analysis_layer="pca"
iteration_num="0"
balance_bool="No"
n_dm="0"
mode_embedding="PCA"
layer_embedding="X_pca"

# Dataset specific parameters
data_dir="${root}/data/synthetic/cluster"
data_file="${data_dir}/cluster_PCA_0.h5ad"
pops="M1 M2 M3"
batch_vec="0 0.75 1 1.25 1.5"
k=30
resolution=0.2
beta=33
downsample=3
pop_col="celltype"

# R methods to run
methods="milo daseq cydar louvain"

# Set default values for environment variables if not set
if [ -z "${SEEDS+x}" ]; then
    SEEDS=(43 44 45)
fi

if [ -z "${ENRICHMENT_VALUES+x}" ]; then
    ENRICHMENT_VALUES=(0.75 0.85 0.95)
fi

# Job counter for array jobs
job_number=0

# Loop through parameters
for p in $pops; do
    for seed in ${SEEDS[@]}; do
        for enr in ${ENRICHMENT_VALUES[@]}; do
            for batch_sd_num in $batch_vec; do
                for method in $methods; do
                    for iteration in $(seq 0 1 $iteration_num); do
                        ((job_number++))
                        
                        # Skip if not the current array task
                        if [ -z "$SLURM_ARRAY_TASK_ID" ] || [ "$job_number" -ne "$SLURM_ARRAY_TASK_ID" ]; then
                            continue
                        fi
                        
                        # Create job ID and paths
                        # Use unified jobid for label input (embedding-independent)
                        unified_jobid=${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}
                        # Use embedding-specific jobid for benchmark output
                        jobid=${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}
                        input_path=${data_dir}/${unified_jobid}
                        save_path=${root}/benchmark_${mode_embedding,,}/$([[ "$data_id" =~ ^(linear|branch|cluster)$ ]] && echo "synthetic" || echo "real")/$data_id/${jobid}
                        save_path_iteration=${save_path}/iteration_${iteration}
                        mkdir -p "$save_path_iteration"
                        
                        echo "Running $jobid method=$method..."

                        # R environment is handled by conda/mamba activation
                        
                        # Execute R method
                        Rscript ${{root}}/scripts/run_DA.r \
                            --file_path ${{data_file}} \
                            --pop ${{p}} \
                            --pop_enr ${{enr}} \
                            --pop_column ${{pop_col}} \
                            --ds_type ${{data_id}} \
                            --batch_sd ${{batch_sd_num}} \
                            --input_file ${{input_path}}/ \
                            --package ${{method}} \
                            --seed ${{seed}} \
                            --layer_embedding ${{layer_embedding}} \
                            --k ${{k}} \
                            --resolution ${{resolution}} \
                            --n_dm ${{n_dm}} \
                            --output_dir ${{save_path}}/
                        exit 0
                    done
                done
            done
        done
    done
done

# Define method-specific command functions

# Submit a slurm array job if not already running in slurm
if [ -z "$SLURM_JOB_ID" ]; then
    jobid="benchmark_${data_id}_${mode_embedding,,}"
    cmd="sbatch -J '$jobid' --time=$time --partition=$partition \
    --mem 8g --out '$root/SlurmLog/${jobid}_%N_%A_%a.out' --array=1-$job_number \
    '$script_path'"
    echo "$cmd"
    eval "$cmd"
fi
