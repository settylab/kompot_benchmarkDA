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
from lib.logger import get_logger

# Pipeline constants (same across all datasets)
PIPELINE_CONSTANTS = {
    "n_conditions": 2,
    "n_replicates": 3,
    "n_batches": 2,
    "m": 2,
    "a_logit": 0.5,
    "layer_embedding": "X_pca",
}


def main():
    parser = argparse.ArgumentParser(
        description="Generate labels for benchmarking - uses config defaults with optional filtering"
    )
    parser.add_argument("--dataset", type=str, required=True, help="Dataset name")
    parser.add_argument(
        "--populations", type=str, help="Filter: comma-separated populations"
    )
    parser.add_argument("--seeds", type=str, help="Filter: comma-separated seeds")
    parser.add_argument(
        "--enrichments", type=str, help="Filter: comma-separated enrichments"
    )
    parser.add_argument(
        "--batch-sds", type=str, help="Filter: comma-separated batch SDs"
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip combinations that already have labels",
    )

    args = parser.parse_args()
    logger = get_logger("generate_labels", verbose=True)

    dataset = args.dataset
    if dataset not in DATASET_CONFIGS:
        logger.error(f"Dataset '{dataset}' not found in configuration")
        logger.info(f"Available datasets: {', '.join(DATASET_CONFIGS.keys())}")
        sys.exit(1)

    # Get dataset configuration
    dataset_config = DATASET_CONFIGS[dataset]

    # Use config defaults, apply filters only if specified
    populations = (
        args.populations.split(",") if args.populations else dataset_config["pops"]
    )
    seeds = [int(s) for s in args.seeds.split(",")] if args.seeds else SEEDS
    enrichments = (
        [float(e) for e in args.enrichments.split(",")]
        if args.enrichments
        else ENRICHMENT_VALUES
    )
    batch_sds = (
        [float(b) for b in args.batch_sds.split(",")]
        if args.batch_sds
        else dataset_config["batch_vec"]
    )

    # Get all dataset-specific settings from config
    pop_col = dataset_config.get("pop_col", "celltype")
    n_dm = dataset_config.get("n_dm", 0)

    # Determine dataset type and data file
    if dataset in ["linear", "branch", "cluster"]:
        data_type = "synthetic"
    else:
        data_type = "real"

    data_file = project_root / "data" / data_type / dataset / "{}.h5ad".format(dataset)

    if not data_file.exists():
        logger.error(f"Data file not found: {data_file}")
        sys.exit(1)

    total_combinations = (
        len(populations) * len(seeds) * len(enrichments) * len(batch_sds)
    )

    logger.section(f"Label Generation for '{dataset}'")
    logger.info("Configuration from dataset_config.py:")
    logger.info(f"  Populations: {populations}")
    logger.info(f"  Seeds: {seeds}")
    logger.info(f"  Enrichments: {enrichments}")
    logger.info(f"  Batch SDs: {batch_sds}")
    logger.info(f"  Population column: {pop_col}")
    logger.info(f"  DM components: {n_dm}")
    logger.info(f"Total combinations: {total_combinations}")
    print("")

    success_count = 0
    fail_count = 0
    skip_count = 0
    current_task = 0

    for pop in populations:
        for seed in seeds:
            for enrichment in enrichments:
                for batch_sd in batch_sds:
                    current_task += 1
                    combo_id = "{}-{}-{}-{}-{}".format(
                        dataset, pop, enrichment, seed, batch_sd
                    )

                    # Create output directory
                    data_dir = project_root / "data" / data_type / dataset
                    output_dir = data_dir / "{}-{}-{}-{}-{}".format(
                        dataset, pop, enrichment, seed, batch_sd
                    )

                    # Check if labels already exist
                    if args.skip_existing:
                        coldata_file = (
                            output_dir
                            / "benchmark_{}_pop_{}_enr{}_seed{}.coldata.csv".format(
                                dataset, pop, enrichment, seed
                            )
                        )
                        if coldata_file.exists():
                            skip_count += 1
                            logger.progress(current_task, total_combinations, combo_id, status="SKIP")
                            continue

                    output_dir.mkdir(parents=True, exist_ok=True)
                    logger.progress(current_task, total_combinations, f"Generating labels for {combo_id}...")

                    # Build command using config and pipeline constants
                    cmd = [
                        "python",
                        str(project_root / "bin" / "generate_bm_data.py"),
                        "--file_path",
                        str(data_file),
                        "--pop",
                        pop,
                        "--pop_enr",
                        str(enrichment),
                        "--pop_column",
                        pop_col,
                        "--ds_type",
                        dataset,
                        "--batch_sd",
                        str(batch_sd),
                        "--seed",
                        str(seed),
                        "--n_conditions",
                        str(PIPELINE_CONSTANTS["n_conditions"]),
                        "--n_replicates",
                        str(PIPELINE_CONSTANTS["n_replicates"]),
                        "--n_batches",
                        str(PIPELINE_CONSTANTS["n_batches"]),
                        "--m",
                        str(PIPELINE_CONSTANTS["m"]),
                        "--a_logit",
                        str(PIPELINE_CONSTANTS["a_logit"]),
                        "--layer_embedding",
                        PIPELINE_CONSTANTS["layer_embedding"],
                        "--n_dm",
                        str(n_dm),
                        "--output_dir",
                        str(output_dir) + "/",
                    ]

                    try:
                        result = subprocess.run(
                            cmd,
                            check=True,
                            capture_output=True,
                            text=True,
                            timeout=3600,  # 1 hour timeout for large datasets with DM computation
                        )
                        logger.progress(current_task, total_combinations, combo_id, status="OK")
                        success_count += 1
                    except subprocess.CalledProcessError as e:
                        logger.progress(current_task, total_combinations, combo_id, status="FAIL")
                        if e.stderr:
                            # Show first and last lines of error
                            error_lines = e.stderr.strip().split("\n")
                            if len(error_lines) > 5:
                                logger.error(f"  First error: {error_lines[0]}")
                                logger.error(f"  Last error: {error_lines[-1]}")
                            else:
                                logger.error(f"  Error: {e.stderr[:200]}")
                        fail_count += 1
                    except subprocess.TimeoutExpired:
                        logger.progress(current_task, total_combinations, combo_id, status="TIMEOUT")
                        fail_count += 1

    print("")
    logger.section("Summary")
    if args.skip_existing and skip_count > 0:
        logger.info(
            f"Summary: {success_count} succeeded, {fail_count} failed, {skip_count} skipped (already exist)"
        )
    else:
        logger.info(f"Summary: {success_count} succeeded, {fail_count} failed")

    return 0 if fail_count == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
