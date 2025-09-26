#!/usr/bin/env python3
"""
Status reporting system for BenchmarkDA pipeline.
Tracks completion status of preprocessing, labels, and benchmark runs.
"""

import os
import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# Add project root to Python path
sys.path.append(str(Path(__file__).resolve().parent.parent))

# Import configurations
from config.dataset_config import DATASET_CONFIGS, SEEDS, ENRICHMENT_VALUES
from config.method_config import PYTHON_METHODS, R_METHODS

class BenchmarkStatus:
    def __init__(self, root_dir="."):
        self.root_dir = Path(root_dir)
        self.status_file = self.root_dir / "status.json"

    def scan_filesystem(self):
        """Scan filesystem to determine current status of all components."""
        status = {
            "last_updated": datetime.now().isoformat(),
            "datasets": {},
            "summary": {
                "total_datasets": 0,
                "preprocessing_complete": 0,
                "labels_complete": 0,
                "benchmarks_complete": 0
            }
        }

        # Get all available datasets
        all_datasets = list(DATASET_CONFIGS.keys())

        for dataset in all_datasets:
            dataset_status = self._check_dataset_status(dataset)
            status["datasets"][dataset] = dataset_status

            # Update summary counts
            status["summary"]["total_datasets"] += 1
            if dataset_status["preprocessing"]["dm"]["complete"] and dataset_status["preprocessing"]["pca"]["complete"]:
                status["summary"]["preprocessing_complete"] += 1
            if dataset_status["labels"]["complete"]:
                status["summary"]["labels_complete"] += 1
            if dataset_status["benchmarks"]["python"]["complete"] and dataset_status["benchmarks"]["r"]["complete"]:
                status["summary"]["benchmarks_complete"] += 1

        return status

    def _check_dataset_status(self, dataset):
        """Check status for a specific dataset."""
        config = DATASET_CONFIGS[dataset]

        # Determine data paths
        if dataset in ['linear', 'branch', 'cluster']:
            data_dir = self.root_dir / "data" / "synthetic" / dataset
            result_type = "synthetic"
        else:
            data_dir = self.root_dir / "data" / "real" / dataset
            result_type = "real"

        dataset_status = {
            "type": result_type,
            "config": {
                "n_dm": self._get_n_dm_for_dataset(dataset),
                "pops": config["pops"],
                "seeds": SEEDS,
                "enrichment_values": ENRICHMENT_VALUES,
                "batch_vec": config["batch_vec"]
            },
            "data_file": {
                "exists": (data_dir / f"{dataset}.h5ad").exists(),
                "path": str(data_dir / f"{dataset}.h5ad")
            },
            "preprocessing": self._check_preprocessing_status(dataset, data_dir),
            "labels": self._check_labels_status(dataset, data_dir, config),
            "benchmarks": self._check_benchmarks_status(dataset, result_type, config)
        }

        return dataset_status

    def _check_preprocessing_status(self, dataset, data_dir):
        """Check preprocessing status for embeddings."""
        n_dm = self._get_n_dm_for_dataset(dataset)

        preprocessing_status = {
            "dm": {
                "required": n_dm > 0,
                "complete": False,
                "files": []
            },
            "pca": {
                "required": True,
                "complete": False,
                "files": []
            }
        }

        # Check for DM files (if required)
        if n_dm > 0:
            dm_file = data_dir / f"{dataset}_DM_{n_dm}.h5ad"
            preprocessing_status["dm"]["complete"] = dm_file.exists()
            if dm_file.exists():
                preprocessing_status["dm"]["files"].append(str(dm_file))
        else:
            preprocessing_status["dm"]["complete"] = True  # Not required

        # Check for PCA files
        pca_file = data_dir / f"{dataset}_PCA_0.h5ad"
        preprocessing_status["pca"]["complete"] = pca_file.exists()
        if pca_file.exists():
            preprocessing_status["pca"]["files"].append(str(pca_file))

        return preprocessing_status

    def _check_labels_status(self, dataset, data_dir, config):
        """Check label generation status."""
        labels_status = {
            "complete": False,
            "total_combinations": 0,
            "completed_combinations": 0,
            "missing_combinations": [],
            "files": []
        }

        # Calculate total combinations needed
        total_combinations = (
            len(config["pops"]) *
            len(SEEDS) *
            len(ENRICHMENT_VALUES) *
            len(config["batch_vec"])
        )
        labels_status["total_combinations"] = total_combinations

        completed_count = 0
        for pop in config["pops"]:
            for seed in SEEDS:
                for enr in ENRICHMENT_VALUES:
                    for batch_sd in config["batch_vec"]:
                        # Unified label naming (embedding-independent)
                        label_id = f"{dataset}-{pop}-{enr}-{seed}-{batch_sd}-No"
                        label_dir = data_dir / label_id

                        if label_dir.exists() and any(label_dir.iterdir()):
                            completed_count += 1
                            labels_status["files"].append(str(label_dir))
                        else:
                            labels_status["missing_combinations"].append(label_id)

        labels_status["completed_combinations"] = completed_count
        labels_status["complete"] = total_combinations > 0 and completed_count == total_combinations

        return labels_status

    def _check_benchmarks_status(self, dataset, result_type, config):
        """Check benchmark completion status."""
        benchmark_status = {
            "python": self._check_method_type_status(dataset, result_type, config, "python"),
            "r": self._check_method_type_status(dataset, result_type, config, "r")
        }

        return benchmark_status

    def _check_method_type_status(self, dataset, result_type, config, method_type):
        """Check status for a specific method type (python or r)."""
        if method_type == "python":
            methods = PYTHON_METHODS
        else:
            methods = R_METHODS

        benchmark_dir = self.root_dir / "benchmark" / result_type / dataset

        method_status = {
            "complete": False,
            "methods": {},
            "total_jobs": 0,
            "completed_jobs": 0,
            "embeddings": ["dm", "pca"]
        }

        total_jobs = 0
        completed_jobs = 0

        for method_name in methods:
            method_info = {
                "complete": False,
                "embeddings": {
                    "dm": {"complete": False, "files": []},
                    "pca": {"complete": False, "files": []}
                }
            }

            for embedding in ["dm", "pca"]:
                embedding_jobs = 0
                embedding_completed = 0

                for pop in config["pops"]:
                    for seed in SEEDS:
                        for enr in ENRICHMENT_VALUES:
                            for batch_sd in config["batch_vec"]:
                                # Job ID with embedding
                                job_id = f"{dataset}-{pop}-{enr}-{seed}-{batch_sd}-No-{embedding}"
                                job_dir = benchmark_dir / job_id / "iteration_0"

                                embedding_jobs += 1
                                total_jobs += 1

                                # Check if results exist
                                if job_dir.exists() and any(job_dir.iterdir()):
                                    # Look for method-specific result files
                                    method_results = list(job_dir.glob(f"*{method_name}*"))
                                    if method_results:
                                        embedding_completed += 1
                                        completed_jobs += 1
                                        method_info["embeddings"][embedding]["files"].extend([str(f) for f in method_results])

                # Only consider complete if there were jobs to do AND they're all done
                method_info["embeddings"][embedding]["complete"] = (
                    embedding_jobs > 0 and embedding_completed == embedding_jobs
                )

            # Method is complete if both embeddings are complete
            method_info["complete"] = (
                method_info["embeddings"]["dm"]["complete"] and
                method_info["embeddings"]["pca"]["complete"]
            )
            method_status["methods"][method_name] = method_info

        method_status["total_jobs"] = total_jobs
        method_status["completed_jobs"] = completed_jobs
        method_status["complete"] = total_jobs > 0 and completed_jobs == total_jobs

        return method_status

    def _get_n_dm_for_dataset(self, dataset):
        """Get number of DM components for dataset."""
        # Try to get from config, fallback to hardcoded values
        try:
            from config.dataset_config import get_n_dm_for_dataset
            return get_n_dm_for_dataset(dataset)
        except ImportError:
            # Fallback hardcoded values
            if dataset in ['linear', 'branch', 'cluster']:
                return 10
            elif dataset in ['covid19-pbmc', 'pancreas']:
                return 30
            elif dataset in ['bcr-xl', 'levine32']:
                return 5
            else:
                return 10

    def save_status(self, status):
        """Save status to JSON file."""
        with open(self.status_file, 'w') as f:
            json.dump(status, f, indent=2)

    def load_status(self):
        """Load status from JSON file."""
        if self.status_file.exists():
            with open(self.status_file, 'r') as f:
                return json.load(f)
        return None

    def print_summary(self, status):
        """Print a summary status report."""
        print("🔍 BenchmarkDA Status Report")
        print("=" * 50)
        print(f"Last Updated: {status['last_updated']}")
        print(f"Total Datasets: {status['summary']['total_datasets']}")
        print(f"Preprocessing Complete: {status['summary']['preprocessing_complete']}/{status['summary']['total_datasets']}")
        print(f"Labels Complete: {status['summary']['labels_complete']}/{status['summary']['total_datasets']}")
        print(f"All Benchmarks Complete: {status['summary']['benchmarks_complete']}/{status['summary']['total_datasets']}")
        print("  (Note: Individual method progress shown in detailed view)")
        print()

    def print_detailed_status(self, status):
        """Print detailed status for each dataset."""
        for dataset, info in status["datasets"].items():
            print(f"📊 Dataset: {dataset} ({info['type']})")
            print("-" * 40)

            # Data file status
            data_icon = "✅" if info["data_file"]["exists"] else "❌"
            print(f"  Data File: {data_icon} {info['data_file']['path']}")

            # Preprocessing status
            dm_icon = "✅" if info["preprocessing"]["dm"]["complete"] else "❌" if info["preprocessing"]["dm"]["required"] else "➖"
            pca_icon = "✅" if info["preprocessing"]["pca"]["complete"] else "❌"
            print(f"  Preprocessing: DM {dm_icon} | PCA {pca_icon}")

            # Labels status
            labels_icon = "✅" if info["labels"]["complete"] else "❌"
            labels_pct = (info["labels"]["completed_combinations"] / info["labels"]["total_combinations"] * 100) if info["labels"]["total_combinations"] > 0 else 0
            print(f"  Labels: {labels_icon} {info['labels']['completed_combinations']}/{info['labels']['total_combinations']} ({labels_pct:.1f}%)")

            # Benchmarks status - show individual methods
            print("  Benchmarks:")

            # Python methods
            python_methods = info["benchmarks"]["python"]["methods"]
            print("    Python:")
            for method_name, method_info in python_methods.items():
                method_icon = "✅" if method_info["complete"] else "❌"
                dm_complete = method_info["embeddings"]["dm"]["complete"]
                pca_complete = method_info["embeddings"]["pca"]["complete"]
                dm_icon = "✅" if dm_complete else "❌"
                pca_icon = "✅" if pca_complete else "❌"
                print(f"      {method_name}: {method_icon} (DM {dm_icon} | PCA {pca_icon})")

            # R methods
            r_methods = info["benchmarks"]["r"]["methods"]
            print("    R:")
            for method_name, method_info in r_methods.items():
                method_icon = "✅" if method_info["complete"] else "❌"
                dm_complete = method_info["embeddings"]["dm"]["complete"]
                pca_complete = method_info["embeddings"]["pca"]["complete"]
                dm_icon = "✅" if dm_complete else "❌"
                pca_icon = "✅" if pca_complete else "❌"
                print(f"      {method_name}: {method_icon} (DM {dm_icon} | PCA {pca_icon})")

            print()

    def print_missing_items(self, status):
        """Print what still needs to be run."""
        print("🚧 Missing Items:")
        print("=" * 50)

        for dataset, info in status["datasets"].items():
            missing_items = []

            if not info["data_file"]["exists"]:
                missing_items.append("Data file")

            if info["preprocessing"]["dm"]["required"] and not info["preprocessing"]["dm"]["complete"]:
                missing_items.append("DM preprocessing")

            if not info["preprocessing"]["pca"]["complete"]:
                missing_items.append("PCA preprocessing")

            if not info["labels"]["complete"]:
                missing_items.append(f"Labels ({info['labels']['total_combinations'] - info['labels']['completed_combinations']} combinations)")

            if not info["benchmarks"]["python"]["complete"]:
                incomplete_methods = [m for m, status in info["benchmarks"]["python"]["methods"].items() if not status["complete"]]
                if incomplete_methods:
                    missing_items.append(f"Python methods: {', '.join(incomplete_methods)}")

            if not info["benchmarks"]["r"]["complete"]:
                incomplete_methods = [m for m, status in info["benchmarks"]["r"]["methods"].items() if not status["complete"]]
                if incomplete_methods:
                    missing_items.append(f"R methods: {', '.join(incomplete_methods)}")

            if missing_items:
                print(f"  {dataset}: {', '.join(missing_items)}")

        print()

    def suggest_commands(self, status):
        """Suggest CLI commands to complete missing work."""
        print("💡 Suggested Commands:")
        print("=" * 50)

        # Find datasets that need preprocessing
        datasets_need_preprocessing = []
        datasets_need_labels = []
        datasets_need_python_benchmarks = []
        datasets_need_r_benchmarks = []

        for dataset, info in status["datasets"].items():
            if not info["data_file"]["exists"]:
                print(f"# Download {dataset} dataset first")
                continue

            needs_preprocessing = (
                (info["preprocessing"]["dm"]["required"] and not info["preprocessing"]["dm"]["complete"]) or
                not info["preprocessing"]["pca"]["complete"]
            )
            if needs_preprocessing:
                datasets_need_preprocessing.append(dataset)

            if not info["labels"]["complete"]:
                datasets_need_labels.append(dataset)

            if not info["benchmarks"]["python"]["complete"]:
                datasets_need_python_benchmarks.append(dataset)

            if not info["benchmarks"]["r"]["complete"]:
                datasets_need_r_benchmarks.append(dataset)

        if datasets_need_preprocessing:
            datasets_str = ",".join(datasets_need_preprocessing)
            print(f"./cli.sh --datasets {datasets_str} preprocess")

        if datasets_need_labels:
            datasets_str = ",".join(datasets_need_labels)
            print(f"./cli.sh --datasets {datasets_str} labels")

        if datasets_need_python_benchmarks:
            datasets_str = ",".join(datasets_need_python_benchmarks)
            print(f"./cli.sh --datasets {datasets_str} --methods python benchmark")

        if datasets_need_r_benchmarks:
            datasets_str = ",".join(datasets_need_r_benchmarks)
            print(f"./cli.sh --datasets {datasets_str} --methods r benchmark")

        print()

def main():
    parser = argparse.ArgumentParser(description="BenchmarkDA Status Report")
    parser.add_argument("--detailed", action="store_true", help="Show detailed status for each dataset")
    parser.add_argument("--missing", action="store_true", help="Show missing items only")
    parser.add_argument("--commands", action="store_true", help="Suggest CLI commands to complete missing work")
    parser.add_argument("--save", action="store_true", help="Save status to JSON file")
    parser.add_argument("--dataset", type=str, help="Show status for specific dataset only")

    args = parser.parse_args()

    # Create status reporter
    reporter = BenchmarkStatus()

    # Scan current status
    status = reporter.scan_filesystem()

    # Filter by dataset if specified
    if args.dataset:
        if args.dataset in status["datasets"]:
            filtered_status = {
                "last_updated": status["last_updated"],
                "datasets": {args.dataset: status["datasets"][args.dataset]},
                "summary": {"total_datasets": 1, "preprocessing_complete": 0, "labels_complete": 0, "benchmarks_complete": 0}
            }
            # Recalculate summary for single dataset
            info = filtered_status["datasets"][args.dataset]
            if info["preprocessing"]["dm"]["complete"] and info["preprocessing"]["pca"]["complete"]:
                filtered_status["summary"]["preprocessing_complete"] = 1
            if info["labels"]["complete"]:
                filtered_status["summary"]["labels_complete"] = 1
            if info["benchmarks"]["python"]["complete"] and info["benchmarks"]["r"]["complete"]:
                filtered_status["summary"]["benchmarks_complete"] = 1
            status = filtered_status
        else:
            print(f"Dataset '{args.dataset}' not found.")
            return

    # Save status if requested
    if args.save:
        reporter.save_status(status)
        print(f"Status saved to {reporter.status_file}")
        print()

    # Print reports based on arguments
    if args.missing:
        reporter.print_missing_items(status)
    elif args.commands:
        reporter.suggest_commands(status)
    elif args.detailed:
        reporter.print_summary(status)
        reporter.print_detailed_status(status)
    else:
        reporter.print_summary(status)
        reporter.print_missing_items(status)
        reporter.suggest_commands(status)

if __name__ == "__main__":
    main()