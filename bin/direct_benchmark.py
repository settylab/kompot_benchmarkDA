#!/usr/bin/env python3
"""
Direct benchmark execution script for the benchmarkDA project.
Eliminates the need for script generation by directly calling methods based on configuration.
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

# Add project root to Python path
sys.path.append(str(Path(__file__).resolve().parent.parent))

# Import configurations
from config.dataset_config import DATASET_CONFIGS, SEEDS, ENRICHMENT_VALUES
from config.method_config import PYTHON_METHODS, R_METHODS

def run_python_method(method_name, method_config, dataset_config, params):
    """Execute a Python method directly."""
    script_path = Path("python_method") / method_config["script"]

    cmd = [
        "python", str(script_path),
        "--file_path", params["data_file"],
        "--pop", params["pop"],
        "--pop_enr", str(params["enrichment"]),
        "--pop_column", dataset_config["pop_col"],
        "--ds_type", params["dataset"],
        "--batch_sd", str(params["batch_sd"]),
        "--input_file", f"{params['input_path']}/",
        "--package", method_name,
        "--seed", str(params["seed"]),
        "--layer_embedding", params["layer_embedding"],
        "--output_dir", f"{params['output_path']}/"
    ]

    # Add method-specific parameters
    if method_config["script"] == "Mellon_bm.py":
        method_params = method_config["params"]
        cmd.extend([
            "--n_dm", str(params["n_dm"]),
            "--mellon_d_method", method_params.get("mellon_d_method", "fractal"),
            "--norm_density", method_params.get("norm_density", "No"),
            "--hyperparameter", method_params.get("hyperparameter", "Yes"),
            "--corrected", method_params.get("corrected", "No"),
            "--ls_factor", str(method_params.get("ls_factor", 1.5)),
            "--ls_mode", params["mode_embedding"]
        ])
    elif method_config["script"] == "meld_bm.py":
        method_params = method_config["params"]
        beta_value = method_params.get("beta", dataset_config["beta"])
        cmd.extend([
            "--beta", str(beta_value),
            "--k_meld", str(dataset_config["k"])
        ])
    elif method_config["script"] == "CNA_bm.py":
        cmd.extend([
            "--k_cna", str(dataset_config["k"])
        ])

    print(f"Executing: {' '.join(cmd)}")
    return subprocess.run(cmd, cwd="python_method")

def run_r_method(method_name, dataset_config, params):
    """Execute an R method directly."""
    cmd = [
        "Rscript", "scripts/run_DA.r",
        "--file_path", params["data_file"],
        "--pop", params["pop"],
        "--pop_enr", str(params["enrichment"]),
        "--pop_column", dataset_config["pop_col"],
        "--ds_type", params["dataset"],
        "--batch_sd", str(params["batch_sd"]),
        "--input_file", f"{params['input_path']}/",
        "--package", method_name,
        "--seed", str(params["seed"]),
        "--layer_embedding", params["layer_embedding"],
        "--k", str(dataset_config["k"]),
        "--resolution", str(dataset_config["resolution"]),
        "--n_dm", str(params["n_dm"]),
        "--output_dir", f"{params['output_path']}/"
    ]

    print(f"Executing: {' '.join(cmd)}")
    return subprocess.run(cmd)

def main():
    """Main function to execute benchmarks directly."""
    parser = argparse.ArgumentParser(description="Direct benchmark execution")

    # Required arguments
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name")
    parser.add_argument("--method_type", type=str, required=True, choices=["python", "r"], help="Method type")
    parser.add_argument("--embeddings", type=str, required=True, help="Comma-separated list of embeddings")
    parser.add_argument("--n_dm", type=int, required=True, help="Number of diffusion map components")

    # Optional arguments
    parser.add_argument("--iteration_num", type=int, default=0, help="Number of iterations")
    parser.add_argument("--balance", type=str, default="No", help="Balance flag")
    parser.add_argument("--methods", type=str, nargs="+", help="Specific methods to run")
    parser.add_argument("--slurm", action="store_true", help="Submit as SLURM array job")

    args = parser.parse_args()

    # Get dataset configuration
    dataset_config = DATASET_CONFIGS.get(args.dataset)
    if not dataset_config:
        raise ValueError(f"Dataset {args.dataset} not found in configuration")

    # Parse embeddings
    embeddings = [e.strip() for e in args.embeddings.split(',')]

    # Get methods to run
    if args.method_type == "python":
        available_methods = PYTHON_METHODS
        if args.methods:
            methods_to_run = {k: v for k, v in PYTHON_METHODS.items() if k in args.methods}
        else:
            methods_to_run = PYTHON_METHODS
    else:
        available_methods = R_METHODS
        if args.methods:
            methods_to_run = {k: v for k, v in R_METHODS.items() if k in args.methods}
        else:
            methods_to_run = R_METHODS

    # Set data paths
    if args.dataset in ['linear', 'branch', 'cluster']:
        data_file = f"data/synthetic/{args.dataset}/{args.dataset}.h5ad"
        data_dir = f"data/synthetic/{args.dataset}"
        data_type = "synthetic"
    else:
        data_file = f"data/real/{args.dataset}/{args.dataset}.h5ad"
        data_dir = f"data/real/{args.dataset}"
        data_type = "real"

    if args.slurm:
        # Generate and submit SLURM job
        generate_slurm_job(args, dataset_config, embeddings, methods_to_run, data_file, data_dir, data_type)
    else:
        # Execute directly
        execute_benchmarks(args, dataset_config, embeddings, methods_to_run, data_file, data_dir, data_type)

def execute_benchmarks(args, dataset_config, embeddings, methods_to_run, data_file, data_dir, data_type):
    """Execute benchmarks directly without SLURM."""

    for embedding in embeddings:
        # Set embedding-specific variables
        if embedding == "dm":
            analysis_layer = "dm"
            layer_embedding = "X_pca"
            mode_embedding = "DM"
        else:
            analysis_layer = "pca"
            layer_embedding = "X_pca"
            mode_embedding = "PCA"

        print(f"Processing embedding: {embedding} ({mode_embedding})")

        # Execute for each parameter combination
        for pop in dataset_config['pops']:
            for seed in SEEDS:
                for enrichment in ENRICHMENT_VALUES:
                    for batch_sd in dataset_config['batch_vec']:
                        for method_name in methods_to_run:
                            for iteration in range(args.iteration_num + 1):
                                # Create paths
                                unified_jobid = f"{args.dataset}-{pop}-{enrichment}-{seed}-{batch_sd}-{args.balance}"
                                embedding_jobid = f"{unified_jobid}-{analysis_layer}"
                                input_path = f"{data_dir}/{unified_jobid}"
                                output_path = f"benchmark/{data_type}/{args.dataset}/{embedding_jobid}/iteration_{iteration}"

                                # Create output directory
                                os.makedirs(output_path, exist_ok=True)

                                # Prepare parameters
                                params = {
                                    "dataset": args.dataset,
                                    "pop": pop,
                                    "enrichment": enrichment,
                                    "seed": seed,
                                    "batch_sd": batch_sd,
                                    "n_dm": args.n_dm,
                                    "data_file": data_file,
                                    "input_path": input_path,
                                    "output_path": output_path,
                                    "layer_embedding": layer_embedding,
                                    "mode_embedding": mode_embedding
                                }

                                print(f"Running {embedding_jobid} method={method_name}...")

                                # Execute method
                                if args.method_type == "python":
                                    method_config = methods_to_run[method_name]
                                    result = run_python_method(method_name, method_config, dataset_config, params)
                                else:
                                    result = run_r_method(method_name, dataset_config, params)

                                if result.returncode != 0:
                                    print(f"Warning: {method_name} failed for {embedding_jobid}")

def generate_slurm_job(args, dataset_config, embeddings, methods_to_run, data_file, data_dir, data_type):
    """Generate a SLURM array job for the benchmarks."""
    # This would generate a SLURM script similar to the old approach but more streamlined
    # For now, just execute directly
    print("SLURM submission not implemented yet, executing directly...")
    execute_benchmarks(args, dataset_config, embeddings, methods_to_run, data_file, data_dir, data_type)

if __name__ == "__main__":
    main()