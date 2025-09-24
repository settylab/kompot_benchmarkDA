#!/bin/bash

# Check if module command is available and load modules on systems that support it
if command -v module &> /dev/null; then
    module purge
    module load ImageMagick/7.1.0-53-GCCcore-12.2.0 || true
    module load GSL/2.7-GCCcore-12.2.0 || true
    module load cuDNN/8.4.1.50-CUDA-11.7.0 || true
fi

# Set up micromamba environment
eval "$(micromamba shell hook --shell bash 2>/dev/null)" || echo "micromamba not available, assuming environment is already activated"

# Set slurm parameters
time=1-00:00:00
partition=campus-new

script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"
if [ -z "${root+x}" ]; then
    export root="$(readlink -f "$script_dir/..")"
fi

cd ${root}/python_method
echo "Working directory: $root/python_method"

# Parameters from command line
data_id="linear"
analysis_layer="pca"
iteration_num="0"
balance_bool="No"
n_dm="0"
mode_embedding="PCA"

# Dataset specific parameters
data_dir="${root}/data/synthetic/linear"
data_file="${data_dir}/linear_PCA_0.h5ad"
layer_embedding="X_pca"
pops="M1 M2 M3 M4 M5 M6 M7"
batch_vec="0 0.75 1 1.25 1.5"
k=30
resolution=1
beta=71
downsample=3
pop_col="celltype"

# Python methods to run
methods="kompot"

# Set default values if environment variables are not set
if [ -z "${SEEDS+x}" ]; then
    SEEDS=("0" "1" "2" "3" "4")
fi

if [ -z "${ENRICHMENT_VALUES+x}" ]; then
    ENRICHMENT_VALUES=("1.5" "2.0" "2.5" "3.0")
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
                        jobid=${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}
                        save_path=${root}/benchmark_${mode_embedding,,}/$([[ "$data_id" =~ ^(linear|branch|cluster)$ ]] && echo "synthetic" || echo "real")/$data_id/${jobid}
                        save_path_iteration=${save_path}/iteration_${iteration}
                        mkdir -p "$save_path_iteration"
                        
                        echo "Running $jobid method=$method..."
                        
                        # Activate environment (with error handling for non-slurm environments)
                        micromamba deactivate 2>/dev/null || true
                        micromamba activate diffabundance 2>/dev/null || echo "Using existing environment"
                        
                        # Execute the appropriate method based on selection
                        if [[ "$method" == "kompot" ]]; then
                            kompot_bm_cmd
                        fi

                        exit 0
                    done
                done
            done
        done
    done
done

# Define method-specific command functions

kompot_bm_cmd() {
    python kompot_bm.py \
        --file_path ${data_file} \
        --pop ${p} \
        --pop_enr ${enr} \
        --pop_column ${pop_col} \
        --ds_type ${data_id} \
        --batch_sd ${batch_sd_num} \
        --input_file ${data_dir}/${jobid}/ \
        --package kompot \
        --seed ${seed} \
        --layer_embedding ${layer_embedding} \
        --n_dm ${n_dm} \
        --ls_factor 10.0 \
        --log_fold_change_threshold 1.0 \
        --ptp_threshold 0.05 \
        --output_dir ${save_path}/
}

# Submit a slurm array job if not already running in slurm
if [ -z "$SLURM_JOB_ID" ]; then
    jobid="benchmark_${data_id}_${mode_embedding,,}"
    cmd="sbatch -J '$jobid' --time=$time --partition=$partition \
    --mem 8g --out '$root/SlurmLog/${jobid}_%N_%A_%a.out' --array=1-$job_number \
    '$script_path'"
    echo "$cmd"
    eval "$cmd"
fi
