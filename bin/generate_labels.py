#!/usr/bin/env python3
"""
Unified label generation utility for BenchmarkDA.
Uses ONLY configured values from dataset_config.py with optional filtering.
Replaces: modified_benchmarkda.sh, modified_benchmarkda_dm_all.sh, generate_labels_filtered.py
"""

import sys
import argparse
import subprocess
import os
from pathlib import Path

# Add project root to path
project_root = Path(__file__).resolve().parent.parent
sys.path.append(str(project_root))

from config.dataset_config import DATASET_CONFIGS, SEEDS, ENRICHMENT_VALUES

# Pipeline constants (same across all datasets)
PIPELINE_CONSTANTS = {
    "n_conditions": 2,
    "n_replicates": 3,
    "n_batches": 2,
    "condition_balance": 1,
    "m": 2,
    "a_logit": 0.5,
    "mode_embedding": "PCA",
    "layer_embedding": "X_pca",
    "balance": "No"
}

def main():
    parser = argparse.ArgumentParser(
        description="Generate labels for benchmarking - uses config defaults with optional filtering"
    )
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name")
    parser.add_argument("--populations", type=str, help="Filter: comma-separated populations")
    parser.add_argument("--seeds", type=str, help="Filter: comma-separated seeds")
    parser.add_argument("--enrichments", type=str, help="Filter: comma-separated enrichments")
    parser.add_argument("--batch-sds", type=str, help="Filter: comma-separated batch SDs")

    args = parser.parse_args()

    dataset = args.dataset
    if dataset not in DATASET_CONFIGS:
        print("Error: Dataset '{}' not found in configuration".format(dataset))
        print("Available datasets: {}".format(", ".join(DATASET_CONFIGS.keys())))
        sys.exit(1)

    # Get dataset configuration
    dataset_config = DATASET_CONFIGS[dataset]

    # Use config defaults, apply filters only if specified
    populations = args.populations.split(',') if args.populations else dataset_config['pops']
    seeds = [int(s) for s in args.seeds.split(',')] if args.seeds else SEEDS
    enrichments = [float(e) for e in args.enrichments.split(',')] if args.enrichments else ENRICHMENT_VALUES
    batch_sds = [float(b) for b in args.batch_sds.split(',')] if args.batch_sds else dataset_config['batch_vec']

    # Get all dataset-specific settings from config
    pop_col = dataset_config.get('pop_col', 'celltype')

    # Determine dataset type and data file
    if dataset in ['linear', 'branch', 'cluster']:
        data_type = 'synthetic'
    else:
        data_type = 'real'

    data_file = project_root / "data" / data_type / dataset / "{}.h5ad".format(dataset)

    if not data_file.exists():
        print("Error: Data file not found: {}".format(data_file))
        sys.exit(1)

    total_combinations = len(populations) * len(seeds) * len(enrichments) * len(batch_sds)

    print("=" * 60)
    print("Label Generation for '{}'".format(dataset))
    print("=" * 60)
    print("Configuration from dataset_config.py:")
    print("  Populations: {}".format(populations))
    print("  Seeds: {}".format(seeds))
    print("  Enrichments: {}".format(enrichments))
    print("  Batch SDs: {}".format(batch_sds))
    print("  Population column: {}".format(pop_col))
    print("")
    print("Total combinations: {}".format(total_combinations))
    print("=" * 60)
    print("")

    success_count = 0
    fail_count = 0

    # Change to python_method directory for script execution
    os.chdir(project_root / "python_method")

    for pop in populations:
        for seed in seeds:
            for enrichment in enrichments:
                for batch_sd in batch_sds:
                    combo_id = "{}-{}-{}-{}-{}".format(dataset, pop, enrichment, seed, batch_sd)

                    # Create output directory
                    data_dir = project_root / "data" / data_type / dataset
                    output_dir = data_dir / "{}-{}-{}-{}-{}-No".format(dataset, pop, enrichment, seed, batch_sd)
                    output_dir.mkdir(parents=True, exist_ok=True)

                    # Build command using config and pipeline constants
                    cmd = [
                        "python", "generate_bm_data.py",
                        "--file_path", str(data_file),
                        "--pop", pop,
                        "--pop_enr", str(enrichment),
                        "--pop_column", pop_col,
                        "--ds_type", dataset,
                        "--batch_sd", str(batch_sd),
                        "--seed", str(seed),
                        "--n_conditions", str(PIPELINE_CONSTANTS["n_conditions"]),
                        "--n_replicates", str(PIPELINE_CONSTANTS["n_replicates"]),
                        "--n_batches", str(PIPELINE_CONSTANTS["n_batches"]),
                        "--condition_balance", str(PIPELINE_CONSTANTS["condition_balance"]),
                        "--m", str(PIPELINE_CONSTANTS["m"]),
                        "--a_logit", str(PIPELINE_CONSTANTS["a_logit"]),
                        "--mode_embedding", PIPELINE_CONSTANTS["mode_embedding"],
                        "--layer_embedding", PIPELINE_CONSTANTS["layer_embedding"],
                        "--balance", PIPELINE_CONSTANTS["balance"],
                        "--output_dir", str(output_dir) + "/"
                    ]

                    try:
                        result = subprocess.run(
                            cmd,
                            check=True,
                            capture_output=True,
                            text=True,
                            timeout=300
                        )
                        print("[OK] {}".format(combo_id))
                        success_count += 1
                    except subprocess.CalledProcessError as e:
                        print("[FAIL] {}".format(combo_id))
                        if e.stderr:
                            # Show first and last lines of error
                            error_lines = e.stderr.strip().split('\n')
                            if len(error_lines) > 5:
                                print("  First error: {}".format(error_lines[0]))
                                print("  Last error: {}".format(error_lines[-1]))
                            else:
                                print("  Error: {}".format(e.stderr[:200]))
                        fail_count += 1
                    except subprocess.TimeoutExpired:
                        print("[TIMEOUT] {}".format(combo_id))
                        fail_count += 1

    print("")
    print("=" * 60)
    print("Summary: {} succeeded, {} failed".format(success_count, fail_count))
    print("=" * 60)

    return 0 if fail_count == 0 else 1

if __name__ == "__main__":
    sys.exit(main())
