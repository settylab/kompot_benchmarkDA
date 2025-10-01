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
    # Get the project root directory
    project_root = Path(__file__).resolve().parent.parent
    script_path = project_root / "python_method" / method_config["script"]

    # Make sure data file path is absolute
    data_file_path = project_root / params["data_file"]

    cmd = [
        "python", str(script_path),
        "--file_path", str(data_file_path),
        "--pop", params["pop"],
        "--pop_enr", str(params["enrichment"]),
        "--pop_column", dataset_config["pop_col"],
        "--ds_type", params["dataset"],
        "--batch_sd", str(params["batch_sd"]),
        "--input_file", f"{str(project_root / params['input_path'])}/",
        "--package", method_name,
        "--seed", str(params["seed"]),
        "--layer_embedding", params["layer_embedding"],
        "--output_dir", f"{str(project_root / params['output_path'])}/"
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
        # CNA methods are disabled due to NumPy 2.0 compatibility issues
        print(f"Warning: CNA method {method_name} is disabled due to compatibility issues")
        return False
    elif method_config["script"] == "kompot_bm.py":
        method_params = method_config["params"]
        cmd.extend([
            "--n_dm", str(params["n_dm"]),
            "--ls_factor", str(method_params.get("ls_factor", 10.0)),
            "--log_fold_change_threshold", str(method_params.get("log_fold_change_threshold", 1.0)),
            "--pvalue_threshold", str(method_params.get("pvalue_threshold", 0.05))
        ])
        if method_params.get("n_landmarks") is not None:
            cmd.extend(["--n_landmarks", str(method_params["n_landmarks"])])
        if method_params.get("force_pca"):
            cmd.append("--force_pca")

    print(f"Executing: {' '.join(cmd)}")
    return subprocess.run(cmd, cwd="python_method", timeout=21600)  # 6 hour timeout

def run_r_method(method_name, dataset_config, params):
    """Execute an R method directly."""
    # Force R methods to use PCA embeddings (n_dm=0) since they don't have DM embeddings
    r_n_dm = 0

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
        "--n_dm", str(r_n_dm),
        "--output_dir", f"{params['output_path']}/"
    ]

    print(f"Executing: {' '.join(cmd)}")

    # Set environment variables for R/renv compatibility in SLURM
    env = os.environ.copy()
    env["RENV_CONFIG_SANDBOX_ENABLED"] = "FALSE"  # Disable renv sandbox to avoid timeouts
    env["R_LIBS_USER"] = ""  # Prioritize conda packages over user library

    # Load GCC module only if needed (for R methods and GLIBC < 2.29)
    needs_gcc_module = False
    if subprocess.run(["which", "module"], capture_output=True).returncode == 0:
        # Check GLIBC version
        try:
            glibc_output = subprocess.run(["ldd", "--version"], capture_output=True, text=True)
            glibc_version = glibc_output.stdout.split('\n')[0].split()[-1]
            major, minor = map(int, glibc_version.split('.'))

            # Only load module if GLIBC < 2.29 (needed for modern R package compilation)
            if major < 2 or (major == 2 and minor < 29):
                needs_gcc_module = True
                print(f"GLIBC {glibc_version} detected - need GCC module for R compilation")
        except:
            # If we can't detect GLIBC, err on the side of caution
            needs_gcc_module = True
            print("Could not detect GLIBC version - loading GCC module as precaution")

    if needs_gcc_module:
        # Check if GCC module is available and load it
        gcc_versions = ["13.3.0", "13.2.0", "12.3.0", "12.2.0", "11.3.0", "11.2.0"]
        for gcc_version in gcc_versions:
            check_cmd = f"module avail GCC/{gcc_version} 2>&1 | grep -q 'GCC/{gcc_version}'"
            if subprocess.run(check_cmd, shell=True, capture_output=True).returncode == 0:
                load_cmd = f"module load GCC/{gcc_version}"
                print(f"Loading GCC/{gcc_version} module for R compilation support")
                cmd = ["/bin/bash", "-c", f"{load_cmd} && {' '.join(cmd)}"]
                break

    return subprocess.run(cmd, timeout=21600, env=env)  # 6 hour timeout

def main():
    """Main function to execute benchmarks directly."""
    parser = argparse.ArgumentParser(description="Direct benchmark execution")

    # Required arguments
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name")
    parser.add_argument("--method_type", type=str, required=True, choices=["python", "r"], help="Method type")
    parser.add_argument("--n_dm", type=int, required=True, help="Number of diffusion map components")
    # Legacy argument for compatibility - embeddings no longer used
    parser.add_argument("--embeddings", type=str, default="pca", help="LEGACY: Ignored")

    # Optional arguments
    parser.add_argument("--balance", type=str, default="No", help="Balance flag")
    parser.add_argument("--methods", type=str, nargs="+", help="Specific methods to run")
    parser.add_argument("--slurm", action="store_true", help="Submit as SLURM array job")
    parser.add_argument("--skip-existing", action="store_true", help="Skip jobs that already have results")

    # Filtering arguments (same as label generation)
    parser.add_argument("--populations", type=str, help="Comma-separated list of populations to run")
    parser.add_argument("--seeds", type=str, help="Comma-separated list of seeds to run")
    parser.add_argument("--enrichments", type=str, help="Comma-separated list of enrichments to run")
    parser.add_argument("--batch-sds", type=str, help="Comma-separated list of batch SDs to run")

    args = parser.parse_args()

    # Get dataset configuration
    dataset_config = DATASET_CONFIGS.get(args.dataset)
    if not dataset_config:
        raise ValueError(f"Dataset {args.dataset} not found in configuration")

    # Apply filtering if specified
    populations = args.populations.split(',') if args.populations else dataset_config['pops']
    seeds = [int(s) for s in args.seeds.split(',')] if args.seeds else SEEDS
    enrichments = [float(e) for e in args.enrichments.split(',')] if args.enrichments else ENRICHMENT_VALUES
    batch_sds = [float(b) for b in getattr(args, 'batch_sds', '').split(',')] if getattr(args, 'batch_sds', None) else dataset_config['batch_vec']

    # Create a filtered config for this run
    filtered_config = dataset_config.copy()
    filtered_config['pops'] = populations
    filtered_config['seeds'] = seeds
    filtered_config['enrichments'] = enrichments
    filtered_config['batch_vec'] = batch_sds

    # all methods use single processing approach
    # embeddings parameter ignored - kept for compatibility

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

    # Get the project root directory
    project_root = Path(__file__).resolve().parent.parent

    # Set data paths
    if args.dataset in ['linear', 'branch', 'cluster']:
        data_file = f"data/synthetic/{args.dataset}/{args.dataset}.h5ad"
        data_dir = str(project_root / f"data/synthetic/{args.dataset}")
        data_type = "synthetic"
    else:
        data_file = f"data/real/{args.dataset}/{args.dataset}.h5ad"
        data_dir = str(project_root / f"data/real/{args.dataset}")
        data_type = "real"

    if args.slurm:
        # Generate and submit SLURM job
        generate_slurm_job(args, filtered_config, methods_to_run, data_file, data_dir, data_type)
    else:
        # Execute directly
        execute_benchmarks(args, filtered_config, methods_to_run, data_file, data_dir, data_type)

def execute_benchmarks(args, dataset_config, methods_to_run, data_file, data_dir, data_type):
    """Execute benchmarks with optional SLURM array support."""

    # Get the project root directory
    project_root = Path(__file__).resolve().parent.parent

    # single processing approach
    analysis_layer = "pca"  # Unified processing layer
    layer_embedding = "X_pca"  # All methods start with PCA embedding
    mode_embedding = "PCA"  # Mode is always PCA-based

    print(f"Processing with PCA-based embeddings")

    # Check if running as SLURM array job
    slurm_task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
    if slurm_task_id:
        print(f"Running as SLURM array task {slurm_task_id}")
        # Execute single task based on SLURM array task ID
        execute_single_task(int(slurm_task_id), args, dataset_config, methods_to_run,
                          data_file, data_dir, data_type, project_root, layer_embedding, mode_embedding)
    else:
        # Execute all parameter combinations locally (using filtered values from config)
        job_number = 0
        for pop in dataset_config['pops']:
            for seed in dataset_config['seeds']:
                for enrichment in dataset_config['enrichments']:
                    for batch_sd in dataset_config['batch_vec']:
                        for method_name in methods_to_run:
                            job_number += 1
                            execute_single_job(job_number, pop, seed, enrichment, batch_sd, method_name,
                                             args, dataset_config, methods_to_run, data_file, data_dir, data_type,
                                             project_root, layer_embedding, mode_embedding)

def execute_single_task(task_id, args, dataset_config, methods_to_run, data_file, data_dir, data_type,
                       project_root, layer_embedding, mode_embedding):
    """Execute a single task for SLURM array job (using filtered config)."""
    # Map task ID to parameter combination
    job_number = 0
    for pop in dataset_config['pops']:
        for seed in dataset_config['seeds']:
            for enrichment in dataset_config['enrichments']:
                for batch_sd in dataset_config['batch_vec']:
                    for method_name in methods_to_run:
                        job_number += 1
                        if job_number == task_id:
                            execute_single_job(job_number, pop, seed, enrichment, batch_sd, method_name,
                                             args, dataset_config, methods_to_run, data_file, data_dir, data_type,
                                             project_root, layer_embedding, mode_embedding)
                            return
    print(f"Error: Task ID {task_id} not found in job range")

def execute_single_job(job_number, pop, seed, enrichment, batch_sd, method_name,
                      args, dataset_config, methods_to_run, data_file, data_dir, data_type,
                      project_root, layer_embedding, mode_embedding):
    """Execute a single benchmark job."""
    # Create job identifiers
    job_id = f"{args.dataset}-{pop}-{enrichment}-{seed}-{batch_sd}-{args.balance}"
    input_path = f"{data_dir}/{job_id}"
    output_path = str(project_root / f"benchmark/{data_type}/{args.dataset}/{job_id}")

    # Check if results already exist
    if args.skip_existing:
        result_file = Path(output_path) / f"benchmark_{args.dataset}_pop_{pop}_enr{enrichment}_seed{seed}_batchEffect{batch_sd}_package_performance.DAresults.{method_name}.csv"
        if result_file.exists():
            print(f"[SKIP] {job_id} method={method_name} (result exists)")
            return

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

    print(f"Running {job_id} method={method_name}...")

    # Execute method
    if args.method_type == "python":
        method_config = methods_to_run[method_name]
        result = run_python_method(method_name, method_config, dataset_config, params)
    else:
        result = run_r_method(method_name, dataset_config, params)

    if result.returncode != 0:
        print(f"Warning: {method_name} failed for {job_id}")

def format_slurm_array(task_ids):
    """Convert a list of task IDs to SLURM array format (e.g., 1,3,5-10,15)."""
    if not task_ids:
        return ""

    task_ids = sorted(task_ids)
    ranges = []
    start = task_ids[0]
    end = task_ids[0]

    for i in range(1, len(task_ids)):
        if task_ids[i] == end + 1:
            # Consecutive, extend range
            end = task_ids[i]
        else:
            # Gap found, save current range
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")
            start = task_ids[i]
            end = task_ids[i]

    # Add final range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{end}")

    return ",".join(ranges)

def generate_slurm_job(args, dataset_config, methods_to_run, data_file, data_dir, data_type):
    """Generate a SLURM array job for the benchmarks."""

    project_root = Path(__file__).resolve().parent.parent

    # Calculate total number of parameter combinations (using filtered config)
    n_pops = len(dataset_config['pops'])
    n_seeds = len(dataset_config['seeds'])
    n_enrichments = len(dataset_config['enrichments'])
    n_batch_vec = len(dataset_config['batch_vec'])
    n_methods = len(methods_to_run)

    total_jobs = n_pops * n_seeds * n_enrichments * n_batch_vec * n_methods

    # If skip-existing is enabled, identify which tasks need to run
    task_list = None
    if args.skip_existing:
        print(f"Detecting incomplete tasks for {args.dataset}...")
        incomplete_tasks = []
        job_number = 0

        for pop in dataset_config['pops']:
            for seed in dataset_config['seeds']:
                for enrichment in dataset_config['enrichments']:
                    for batch_sd in dataset_config['batch_vec']:
                        for method_name in methods_to_run:
                            job_number += 1

                            # Check if result exists
                            job_id = f"{args.dataset}-{pop}-{enrichment}-{seed}-{batch_sd}-{args.balance}"
                            output_path = project_root / f"benchmark/{data_type}/{args.dataset}/{job_id}"
                            result_file = output_path / f"benchmark_{args.dataset}_pop_{pop}_enr{enrichment}_seed{seed}_batchEffect{batch_sd}_package_performance.DAresults.{method_name}.csv"

                            if not result_file.exists():
                                incomplete_tasks.append(job_number)

        if incomplete_tasks:
            task_list = incomplete_tasks
            print(f"  Found {len(incomplete_tasks)} incomplete tasks out of {total_jobs}")
        else:
            print(f"  All {total_jobs} tasks are complete!")
            return  # Nothing to do

    print(f"Generating SLURM array job for {total_jobs} tasks...")
    print(f"  Dataset: {args.dataset}")
    print(f"  Methods: {list(methods_to_run.keys())}")
    print(f"  Populations: {n_pops}")
    print(f"  Seeds: {n_seeds}")
    print(f"  Enrichments: {n_enrichments}")
    print(f"  Batch effects: {n_batch_vec}")

    # Create SLURM script
    script_name = f"slurm_benchmark_{args.dataset}_{args.method_type}.sh"
    script_path = project_root / "benchmark_scripts" / script_name

    # Ensure benchmark_scripts directory exists
    os.makedirs(project_root / "benchmark_scripts", exist_ok=True)

    # Format array specification
    if task_list:
        # Convert list of task IDs to compact SLURM array format
        array_spec = format_slurm_array(task_list)
    else:
        array_spec = f"1-{total_jobs}"

    # Prepare methods list for command line
    methods_arg = ""
    if hasattr(args, 'methods') and args.methods:
        methods_str = " ".join(args.methods)
        methods_arg = f" \\\n        --methods {methods_str}"

    # Don't pass skip-existing to runtime since we already filtered at submission
    skip_arg = ""

    slurm_content = f"""#!/bin/bash
#SBATCH --job-name=bm_{args.dataset}_{args.method_type}
#SBATCH --array={array_spec}
#SBATCH --partition=campus-new
#SBATCH --time=6:00:00
#SBATCH --mem=32G
#SBATCH --output={project_root}/SlurmLog/benchmark_{args.dataset}_{args.method_type}_%A_%a.out
#SBATCH --error={project_root}/SlurmLog/benchmark_{args.dataset}_{args.method_type}_%A_%a.err

# Benchmark execution script
# Generated automatically
# Environment should be pre-activated before submitting this job

set -e

# Environment variables (including module variables) should be inherited from pre-activated environment
# Verify that required environment variables are present
if [ -z "$CONDA_PREFIX" ]; then
    echo "Warning: CONDA_PREFIX not set - mamba environment may not be properly inherited"
fi

# Set additional environment variables for R package management
export R_LIBS_USER=""
export RENV_CONFIG_SANDBOX_ENABLED=FALSE
export MAMBA_NO_BANNER=1
export CONDA_QUIET=1

# Navigate to project root
cd {project_root}

# Execute direct benchmark with SLURM array task ID
# Environment is inherited from CLI via --export=ALL
python bin/direct_benchmark.py \\
    --dataset {args.dataset} \\
    --method_type {args.method_type} \\
    --n_dm {args.n_dm} \\
    --balance {args.balance}{methods_arg}{skip_arg}

echo "Benchmark completed for task $SLURM_ARRAY_TASK_ID"
"""

    # Write SLURM script
    with open(script_path, 'w') as f:
        f.write(slurm_content)

    # Make script executable
    os.chmod(script_path, 0o755)

    print(f"SLURM script created: {script_path}")
    print(f"Submit with: sbatch {script_path}")
    print(f"Monitor with: squeue -u $USER")

    return str(script_path)

if __name__ == "__main__":
    main()