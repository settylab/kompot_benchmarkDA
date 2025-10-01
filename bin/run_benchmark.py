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
from config.dataset_config import DATASET_CONFIGS, SEEDS, ENRICHMENT_VALUES
from config.dataset_config import (
    get_data_file_path,
    get_job_id,
    get_input_dir,
    get_output_dir,
)
from config.method_config import PYTHON_METHODS, R_METHODS
from config.method_config import get_python_method_cmd, get_r_method_cmd


def generate_benchmark_script(args):
    """Generate a standard shell script for running benchmark analyses."""

    # Parse embeddings list
    embeddings = [e.strip() for e in args.embeddings.split(",")]

    # Base script content with user-agnostic environment setup
    script_content = """#!/bin/bash

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

"""

    # Determine working directory based on method type
    if args.method_type == "python":
        script_content += """cd "${root}/python_method"
echo "Working directory: ${root}/python_method"
"""
    else:
        script_content += """cd "${root}/scripts"
echo "Working directory: ${root}/scripts"
"""

    # Get dataset configuration
    dataset_config = DATASET_CONFIGS.get(args.dataset)
    if not dataset_config:
        raise ValueError(f"Dataset {args.dataset} not found in configuration")

    # Add variables from CLI arguments
    script_content += f"""
# Parameters from command line
data_id="{args.dataset}"
iteration_num="{args.iteration_num}"
balance_bool="{args.balance}"
n_dm="{args.n_dm}"
embeddings="{args.embeddings}"

"""

    # Dataset specific variables
    script_content += f"""# Dataset specific parameters
data_dir="${{root}}/data/{('synthetic/' + args.dataset) if args.dataset in ['linear', 'branch', 'cluster'] else args.dataset}"
pops="{' '.join(dataset_config['pops'])}"
batch_vec="{' '.join(map(str, dataset_config['batch_vec']))}"
k={dataset_config['k']}
resolution={dataset_config['resolution']}
beta={dataset_config['beta']}
downsample={dataset_config['downsample']}
pop_col="{dataset_config['pop_col']}"

"""

    # Add methods to use and environment variables
    if args.method_type == "python":
        script_content += f"""# Python methods to run
methods="{' '.join(args.methods) if args.methods else ' '.join(PYTHON_METHODS.keys())}"

# Set default values for environment variables if not set
if [ -z "${{SEEDS+x}}" ]; then
    SEEDS=({' '.join(map(str, SEEDS))})
fi

if [ -z "${{ENRICHMENT_VALUES+x}}" ]; then
    ENRICHMENT_VALUES=({' '.join(map(str, ENRICHMENT_VALUES))})
fi
"""
    else:
        script_content += f"""# R methods to run
methods="{' '.join(args.methods) if args.methods else ' '.join(R_METHODS.keys())}"

# Set default values for environment variables if not set
if [ -z "${{SEEDS+x}}" ]; then
    SEEDS=({' '.join(map(str, SEEDS))})
fi

if [ -z "${{ENRICHMENT_VALUES+x}}" ]; then
    ENRICHMENT_VALUES=({' '.join(map(str, ENRICHMENT_VALUES))})
fi
"""

    # Job counter and nested loops
    script_content += """
# Job counter for array jobs
job_number=0

# Process each embedding type
for embedding in $(echo "$embeddings" | tr ',' ' '); do
    # Set embedding-specific variables
    if [ "$embedding" = "dm" ]; then
        analysis_layer="dm"
        layer_embedding="X_pca"
        mode_embedding="DM"
    else
        analysis_layer="pca"
        layer_embedding="X_pca"
        mode_embedding="PCA"
    fi

    # Set data file path based on dataset type
    if [[ "$data_id" =~ ^(linear|branch|cluster)$ ]]; then
        data_file="${root}/data/synthetic/${data_id}/${data_id}.h5ad"
        data_type="synthetic"
    else
        data_file="${root}/data/real/${data_id}/${data_id}.h5ad"
        data_type="real"
    fi

    echo "Processing embedding: $embedding ($mode_embedding)"

    # Loop through parameters
    for p in $pops; do
        for seed in ${SEEDS[@]}; do
            for enr in ${ENRICHMENT_VALUES[@]}; do
                for batch_sd_num in $batch_vec; do
                    for method in $methods; do
                        for iteration in $(seq 0 1 $iteration_num); do
                            ((job_number++))

                            # Skip if not the current array task
                            if [ -n "$SLURM_ARRAY_TASK_ID" ] && [ "$job_number" -ne "$SLURM_ARRAY_TASK_ID" ]; then
                                continue
                            fi

                            # Create job ID and paths
                            # Use standard jobid for label input (embedding-independent)
                            standard_jobid="${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}"
                            # Use embedding-specific jobid for benchmark output
                            jobid="${data_id}-${p}-${enr}-${seed}-${batch_sd_num}-${balance_bool}-${analysis_layer}"
                            input_path="${data_dir}/${standard_jobid}"
                            save_path="${root}/benchmark/${data_type}/${data_id}/${jobid}"
                            save_path_iteration="${save_path}/iteration_${iteration}"
                            mkdir -p "$save_path_iteration"

                            echo "Running ${jobid} method=${method}..."
"""

    # Method-specific commands
    if args.method_type == "python":
        script_content += """
                            # Activate benchmarkda environment (user-agnostic)
                            activate_benchmarkda_environment || echo "Using existing environment"

                            # Execute the appropriate method based on selection
"""
        # Add Python method conditionals
        for method in PYTHON_METHODS:
            method_config = PYTHON_METHODS[method]
            script_name = method_config["script"].replace(".py", "")
            script_content += f"""                            if [[ "$method" == "{method}" ]]; then
                                {script_name}_cmd
                            fi
"""
        script_content += "                            exit 0\n"
    else:
        script_content += """
                            # R environment is handled by conda/mamba activation
                            activate_benchmarkda_environment || echo "Using existing environment"

                            # Execute R method
                            Rscript "${root}/scripts/run_DA.r" \\
                                --file_path "${data_file}" \\
                                --pop "${p}" \\
                                --pop_enr "${enr}" \\
                                --pop_column "${pop_col}" \\
                                --ds_type "${data_id}" \\
                                --batch_sd "${batch_sd_num}" \\
                                --input_file "${input_path}/" \\
                                --package "${method}" \\
                                --seed "${seed}" \\
                                --layer_embedding "${layer_embedding}" \\
                                --k "${k}" \\
                                --resolution "${resolution}" \\
                                --n_dm "${n_dm}" \\
                                --output_dir "${save_path}/"
                            exit 0
"""

    # Close loops
    script_content += """                        done
                    done
                done
            done
        done
    done
done

# Define method-specific command functions
"""

    # Add method command functions for Python methods
    if args.method_type == "python":
        for method in PYTHON_METHODS:
            method_config = PYTHON_METHODS[method]
            script_name = method_config["script"].replace(".py", "")
            script_content += f"""
{script_name}_cmd() {{
    python "{method_config['script']}" \\
        --file_path "${{data_file}}" \\
        --pop "${{p}}" \\
        --pop_enr "${{enr}}" \\
        --pop_column "${{pop_col}}" \\
        --ds_type "${{data_id}}" \\
        --batch_sd "${{batch_sd_num}}" \\
        --input_file "${{input_path}}/" \\
        --package "{method}" \\
        --seed "${{seed}}" \\
        --layer_embedding "${{layer_embedding}}"""

            # Add method-specific parameters
            if method_config["script"] == "Mellon_bm.py":
                params = method_config["params"]
                script_content += f""" \\
        --n_dm "${{n_dm}}" \\
        --mellon_d_method "{params.get('mellon_d_method', 'fractal')}" \\
        --norm_density "{params.get('norm_density', 'No')}" \\
        --hyperparameter "{params.get('hyperparameter', 'Yes')}" \\
        --corrected "{params.get('corrected', 'No')}" \\
        --ls_factor "{params.get('ls_factor', 1.5)}" \\
        --ls_mode "${{mode_embedding}}"""
            elif method_config["script"] == "meld_bm.py":
                params = method_config["params"]
                beta_value = params.get("beta", "${{beta}}")
                script_content += f""" \\
        --beta "{beta_value}" \\
        --k_meld "${{k}}"""
            elif method_config["script"] == "CNA_bm.py":
                script_content += f""" \\
        --k_cna "${{k}}"""

            script_content += f""" \\
        --output_dir "${{save_path}}/"
}}
"""

    # Slurm job submission
    script_content += """
# Submit a slurm array job if not already running in slurm
if [ -z "$SLURM_JOB_ID" ]; then
    jobid="benchmark_${data_id}"
    cmd="sbatch -J '${jobid}' --time=${time} --partition=${partition} \\
    --mem 8g --out '${root}/SlurmLog/${jobid}_%N_%A_%a.out' --array=1-${job_number} \\
    '${script_path}'"
    echo "$cmd"
    eval "$cmd"
fi
"""

    return script_content


def main():
    """Main function to parse arguments and generate script."""
    parser = argparse.ArgumentParser(description="Generate benchmark scripts")

    # Required arguments
    parser.add_argument(
        "--dataset",
        type=str,
        required=True,
        help="Dataset name (e.g., linear, branch, cluster, covid19-pbmc)",
    )
    parser.add_argument(
        "--method_type",
        type=str,
        required=True,
        choices=["python", "r"],
        help="Method type: python or r",
    )
    parser.add_argument(
        "--embeddings",
        type=str,
        required=True,
        help="Comma-separated list of embeddings (dm,pca)",
    )
    parser.add_argument(
        "--n_dm", type=int, required=True, help="Number of diffusion map components"
    )

    # Optional arguments with defaults
    parser.add_argument(
        "--iteration_num", type=int, default=0, help="Number of iterations (default: 0)"
    )
    parser.add_argument(
        "--balance", type=str, default="No", help="Balance flag (default: No)"
    )
    parser.add_argument(
        "--methods", type=str, nargs="+", help="Methods to run (space-separated list)"
    )
    parser.add_argument("--output", type=str, help="Output script path")

    args = parser.parse_args()

    # Generate script
    script_content = generate_benchmark_script(args)

    # Determine output path
    if args.output:
        output_path = args.output
    else:
        output_path = f"./run_benchmark_{args.dataset}_{args.method_type}.sh"

    # Write script to file
    with open(output_path, "w") as f:
        f.write(script_content)

    # Make executable
    os.chmod(output_path, 0o755)

    print(f"Benchmark script generated: {output_path}")
    print(f"Run with: bash {output_path}")


if __name__ == "__main__":
    main()
