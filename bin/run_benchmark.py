#!/usr/bin/env python3
"""
Unified benchmark runner script for the benchmarkDA project.
Generates shell scripts for running benchmark analyses with consistent parameters.
"""

import os
import sys
import argparse
from pathlib import Path
import json

# Add project root to Python path
sys.path.append(str(Path(__file__).resolve().parent.parent))

# Import configurations
from config.dataset_config import DATASET_CONFIGS_PCA,DATASET_CONFIGS_DM, SEEDS, ENRICHMENT_VALUES
from config.dataset_config import get_data_file_path, get_job_id, get_input_dir, get_output_dir
from config.method_config import PYTHON_METHODS, R_METHODS
from config.method_config import get_python_method_cmd, get_r_method_cmd

def generate_benchmark_script(args):
    """Generate a shell script for running benchmark analyses."""
    
    # Base script content with environment setup
    script_content = """#!/bin/bash


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

"""
    
    # Determine working directory based on method type
    if args.method_type == "python":
        script_content += """cd ${root}/python_method
echo "Working directory: $root/python_method"
"""
    else:
        script_content += """cd ${root}/scripts
echo "Working directory: $root/scripts"
"""

    # Get dataset configuration
    if args.mode_embedding == "PCA":
        dataset_config = DATASET_CONFIGS_PCA.get(args.dataset)
        if not dataset_config:
            raise ValueError(f"Dataset {args.dataset} not found in configuration")
    elif args.mode_embedding == "DM":
        dataset_config = DATASET_CONFIGS_DM.get(args.dataset)
        if not dataset_config:
            raise ValueError(f"Dataset {args.dataset} not found in configuration")
    
    # Add variables from CLI arguments
    script_content += f"""
# Parameters from command line
data_id="{args.dataset}"
analysis_layer="{args.analysis_layer}"
iteration_num={args.iteration_num}
balance_bool="{args.balance}"
n_dm={args.n_dm}
mode_embedding="{args.mode_embedding}"
layer_embedding="{args.layer_embedding}"
batch_sd_mode="{args.batch_sd_mode}"

"""
    
    # Dataset specific variables
    script_content += f"""# Dataset specific parameters
data_dir="${{root}}/data/{('synthetic/' + args.dataset) if args.dataset in ['linear', 'branch', 'cluster'] else ('real/' + args.dataset)}"
data_file_name="${{data_dir}}/{args.dataset}_{args.mode_embedding}_{args.n_dm}"
pops="{' '.join(dataset_config['pops'])}"
batch_vec="{' '.join(map(str, dataset_config['batch_vec_'+args.batch_sd_mode]))}"
k={dataset_config['k']}
resolution={dataset_config['resolution']}
beta={dataset_config['beta']}
downsample={dataset_config['downsample']}
pop_col="{dataset_config['pop_col']}"
seeds="{SEEDS}"
enr_prob="{ENRICHMENT_VALUES}"
"""

    if args.method_type == "python":
        script_content += """
data_file=${data_file_name}.h5ad
"""
    else:
        script_content += """
data_file=${data_file_name}.rds
"""
    # Add methods to use
    if args.method_type == "python":
        script_content += f"""# Python methods to run
methods="{' '.join(args.methods) if args.methods else ' '.join(PYTHON_METHODS.keys())}"
"""
    else:
        script_content += f"""# R methods to run
methods="{' '.join(args.methods) if args.methods else ' '.join(R_METHODS.keys())}"
"""
        
    script_content += """
SEED_OUTPUT=$(python -c "SEED = $seeds; print(' '.join(map(str, SEED)))")
echo "$SEED_OUTPUT"
ENR_OUTPUT=$(python -c "ENRICHMENT_VALUES = $enr_prob; print(' '.join(map(str, ENRICHMENT_VALUES)))")
"""

    script_content += """
# Define method-specific command functions
"""
    # Add method command functions for Python methods
    if args.method_type == "python":
        for method in PYTHON_METHODS:
            method_config = PYTHON_METHODS[method]
            script_content += f"""
{method}_bm_cmd() {{
    python {method_config["script"]} \\
        --file_path ${{data_file}} \\
        --pop ${{p}} \\
        --pop_enr ${{enr}} \\
        --pop_column ${{pop_col}} \\
        --ds_type ${{data_id}} \\
        --batch_sd ${{batch_sd_num}} \\
        --input_file ${{data_dir}}/${{jobid}}/ \\
        --package {method} \\
        --seed ${{seed}} \\
        --layer_embedding ${{layer_embedding}}"""
            
            # Add method-specific parameters
            if method_config["script"] == "Mellon_bm.py":
                params = method_config["params"]
                script_content += f""" \\
        --n_dm 0 \\
        --mellon_d_method "{params.get('mellon_d_method', 'fractal')}" \\
        --norm_density "{params.get('norm_density', 'No')}" \\
        --hyperparameter "{params.get('hyperparameter', 'Yes')}" \\
        --corrected "{params.get('corrected', 'No')}" \\
        --ls_factor {params.get('ls_factor', 1.5)} \\
        --ls_mode ${{mode_embedding}}"""
            elif method_config["script"] == "meld_bm.py":
                params = method_config["params"]
                beta_value = params.get('beta')[args.dataset]
                script_content += f""" \\
        --beta {beta_value} \\
        --k_meld ${{k}}"""
            elif method_config["script"] == "CNA_bm.py":
                script_content += f""" \\
        --k_cna ${{k}}"""
            
            script_content += f""" \\
        --output_dir ${{save_path}}/
}}
"""

    # Job counter and nested loops
    script_content += """
# Job counter for array jobs
job_number=0

# Check whether the dataset exist
if [ -f "$data_file" ]; then
    echo "$data_file exists, start to run benchmarking for the dataset"
# Loop through parameters
    for p in $pops; do
        for seed in ${SEED_OUTPUT[@]}; do
            echo "Using seed: $seed"
            for enr in ${ENR_OUTPUT[@]}; do
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
"""

    # Method-specific commands
    if args.method_type == "python":
        script_content += """                        
                            # Activate environment (with error handling for non-slurm environments)
                            micromamba deactivate 2>/dev/null || true
                            micromamba activate diffabundance 2>/dev/null || echo "Using existing environment"
                            
                            # Execute the appropriate method based on selection
"""
#         script_content += """                        
#                         # Activate environment (with error handling for non-slurm environments)
#                         conda deactivate
#                         conda activate DiffAbundance
                        
#                         # Execute the appropriate method based on selection
# """
        # Add Python method conditional
        method_first = list(PYTHON_METHODS.keys())[0]
        script_content += f"""
                            if [[ "$method" == "{method_first}" ]]; then
                                {method_first}_bm_cmd
                                exit $!
"""
        for method in list(PYTHON_METHODS.keys())[1:-1]:
            method_config = PYTHON_METHODS[method]
            script_content += f"""
                            elif [[ "$method" == "{method}" ]]; then
                                {method}_bm_cmd
                                exit $!
"""
        method_last = list(PYTHON_METHODS.keys())[-1]
        script_content += f"""
                            elif [[ "$method" == "{method_last}" ]]; then
                                {method_last}_bm_cmd
                                exit $!
                            fi
                        done
                    done
                done
            done
        done
    done
"""
        script_content += """
# # Submit a slurm array job if not already running in slurm
    if [ -f ${root}/python_method/prev_job_id_${data_id}_${mode_embedding}.txt ]; then
        prev_job_id=$(cat "${root}/python_method/prev_job_id_${data_id}_${mode_embedding}.txt")
        echo "Child2: Found dependency job ID $prev_job_id"
    fi
    all_statuses=$(sacct -j "$prev_job_id" --format=State --noheader | tr -d ' ')

    # Check if any job info exists for that job id
    if [ -z "$all_statuses" ]; then
        echo "No job found for ID $prev_job_id in Slurm. Submitting without dependency."
        dependency_flag=""
    else
        # Remove any empty lines (if any)
        all_statuses=$(echo "$all_statuses" | grep -v '^$')
        
        # Check for unfinished jobs (i.e. statuses other than COMPLETED, FAILED, or CANCELLED)
        unfinished=$(echo "$all_statuses" | grep -v -E '^(COMPLETED|FAILED|CANCELLED)$')
        
        if [ -z "$unfinished" ]; then
            echo "All dependent jobs ($prev_job_id) are finished. Submitting without dependency."
            dependency_flag=""
        else
            echo "Some dependent jobs are still running or pending."
            dependency_flag="--dependency=afterany:$prev_job_id"
        fi
    fi

    if [ -z "$SLURM_JOB_ID" ]; then
        jobid="benchmark_${data_id}_${mode_embedding,,}"
        cmd="sbatch -J '$jobid' $dependency_flag --time=$time --partition=$partition \\
        --mem 8g --out '$root/SlurmLog/${jobid}_%N_%A_%a.out' --array=1-$job_number \\
        '$script_path'"
        echo "$cmd"
        eval "$cmd"
    fi
"""

        script_content += """
else
    echo "$data_file is not existed, skip for this dataset"
fi
"""
    else:
        script_content += """
                            # Load R environment
                            if command -v module &> /dev/null; then
                                module load R/4.3.1-gfbf-2022b || true
                            fi
                            
                            # Execute R method
"""
        # Add R method code
        script_content += """                            Rscript ${root}/scripts/run_DA.r \\
                                ${data_file} $method $seed $p\\
                                --pop_enrichment $enr \\
                                --data_id ${data_id} \\
                                --batchEffect_sd ${batch_sd_num} \\
                                --data_dir ${data_dir}/${jobid}/iteration_${iteration}/ \\
                                --k ${k} \\
                                --resolution ${resolution} \\
                                --downsample ${downsample} \\
                                --outdir ${save_path_iteration}/
                            exit $!
                        done
                    done
                done
            done
        done
    done
"""
        script_content += """
# # Submit a slurm array job if not already running in slurm
    if [ -f ${root}/python_method/prev_job_id_${data_id}_${mode_embedding}.txt ]; then
        prev_job_id=$(cat "${root}/python_method/prev_job_id_${data_id}_${mode_embedding}.txt")
        echo "Child2: Found dependency job ID $prev_job_id"
    fi

        # Retrieve all job statuses associated with $prev_job_id
    all_statuses=$(sacct -j "$prev_job_id" --format=State --noheader | tr -d ' ')

    # Check if any job info exists for that job id
    if [ -z "$all_statuses" ]; then
        echo "No job found for ID $prev_job_id in Slurm. Submitting without dependency."
        dependency_flag=""
    else
        # Remove any empty lines (if any)
        all_statuses=$(echo "$all_statuses" | grep -v '^$')
        
        # Check for unfinished jobs (i.e. statuses other than COMPLETED, FAILED, or CANCELLED)
        unfinished=$(echo "$all_statuses" | grep -v -E '^(COMPLETED|FAILED|CANCELLED)$')
        
        if [ -z "$unfinished" ]; then
            echo "All dependent jobs ($prev_job_id) are finished. Submitting without dependency."
            dependency_flag=""
        else
            echo "Some dependent jobs are still running or pending."
            dependency_flag="--dependency=afterany:$prev_job_id"
        fi
    fi

    if [ -z "$SLURM_JOB_ID" ]; then
        jobid="benchmark_${data_id}_${mode_embedding,,}"
        cmd="sbatch -J '$jobid' $dependency_flag --time=$time --partition=$partition \\
        --mem 8g --out '$root/SlurmLog/${jobid}_%N_%A_%a.out' --array=1-$job_number \\
        '$script_path'"
        echo "$cmd"
        eval "$cmd"
    fi
    rm "${root}/python_method/prev_job_id_${data_id}_${mode_embedding}.txt"
"""


        script_content +="""
else
    echo "$data_file is not existed, skip for this dataset"
fi
"""

    return script_content

def main():
    """Main function to parse arguments and generate script."""
    parser = argparse.ArgumentParser(description="Generate benchmark scripts")
    
    # Required arguments
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name (e.g., linear, branch, cluster, covid19-pbmc)")
    parser.add_argument("--method_type", type=str, required=True, choices=["python", "r"], help="Method type: python or r")
    
    # Optional arguments with defaults
    parser.add_argument("--analysis_layer", type=str, default="pca", help="Analysis layer name (default: pca)")
    parser.add_argument("--iteration_num", type=int, default=0, help="Number of iterations (default: 0)")
    parser.add_argument("--balance", type=str, default="No", help="Balance flag (default: No)")
    parser.add_argument("--n_dm", type=int, default=0, help="Number of diffusion map components (default: 0)")
    parser.add_argument("--mode_embedding", type=str, default="PCA", choices=["PCA", "DM"], help="Embedding mode (default: PCA)")
    parser.add_argument("--methods", type=str, nargs="+", help="Methods to run (space-separated list)")
    parser.add_argument("--layer_embedding", type=str, default="X_pca", help="the obsm layer for embedding")
    parser.add_argument("--batch_sd_mode", type=str, choices=["modified", "orig"],help="decide whether we want to use the original simulation of batch effect or the modified")
    parser.add_argument("--output", type=str, help="Output script path")
    
    args = parser.parse_args()
    
    # Generate script
    script_content = generate_benchmark_script(args)
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        embedding = args.mode_embedding.lower()
        method_type = args.method_type
        output_path = f"./run_benchmark_{args.dataset}_{embedding}_{method_type}.sh"
    
    # Write script to file
    with open(output_path, "w") as f:
        f.write(script_content)
    
    # Make executable
    os.chmod(output_path, 0o755)
    
    print(f"Benchmark script generated: {output_path}")
    print(f"Run with: bash {output_path}")

if __name__ == "__main__":
    main()