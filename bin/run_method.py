#!/usr/bin/env python3
"""
Unified method runner for benchmarkDA.
Dynamically loads method-specific logic from methods/<method_name>/benchmark.py
This enables easy expandability - just add a new method folder with benchmark.py
"""

import argparse
import sys
from pathlib import Path
import importlib

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib import data_loader, helper_functions
from lib.logger import get_logger


def get_method_module(method_name):
    """
    Dynamically import the benchmark module for a method.

    Args:
        method_name: Name of the method (e.g., 'meld', 'mellon', 'kompot')

    Returns:
        The imported benchmark module

    Raises:
        ImportError: If the method's benchmark module doesn't exist
    """
    # Map method names to their folder names (handle variants like meld_default -> meld)
    base_method = method_name.split("_")[0] if "_" in method_name else method_name

    try:
        module = importlib.import_module(f"methods.{base_method}.benchmark")
        return module
    except ImportError as e:
        raise ImportError(
            f"Could not find benchmark module for method '{method_name}'. "
            f"Expected: methods/{base_method}/benchmark.py"
        ) from e


def main():
    """Main function to run any DA method."""
    # Create base parser with common arguments
    parser = argparse.ArgumentParser(
        description="Run a DA method for benchmarking.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Common arguments (shared across all methods)
    parser.add_argument(
        "--method",
        type=str,
        required=True,
        help="Method name (meld, mellon, kompot, cna, etc.)",
    )
    parser.add_argument(
        "--file_path", type=str, required=True, help="Path to the data file"
    )
    parser.add_argument("--pop", type=str, required=True, help="Population")
    parser.add_argument(
        "--pop_enr", type=float, required=True, help="Population enrichment"
    )
    parser.add_argument(
        "--pop_column", type=str, required=True, help="Population column"
    )
    parser.add_argument("--ds_type", type=str, required=True, help="Type of dataset")
    parser.add_argument(
        "--batch_sd", type=float, required=True, help="Batch standard deviation"
    )
    parser.add_argument(
        "--input_file",
        type=str,
        required=True,
        help="Input file path (label directory)",
    )
    parser.add_argument(
        "--package", type=str, required=True, help="Package name for output files"
    )
    parser.add_argument(
        "--seed", type=int, required=True, help="Seed for random number generation"
    )
    parser.add_argument(
        "--layer_embedding",
        type=str,
        required=True,
        help="Layer embedding (X_pca or DM_EigenVectors)",
    )
    parser.add_argument(
        "--output_dir", type=str, required=True, help="Output directory path"
    )

    # Parse known args first to get the method name
    args, remaining = parser.parse_known_args()
    logger = get_logger("run_method", verbose=True)

    # Load method-specific module and add its arguments
    try:
        method_module = get_method_module(args.method)
        method_module.add_arguments(parser)
    except ImportError as e:
        logger.error(f"{e}")
        sys.exit(1)

    # Parse all arguments (including method-specific ones)
    args = parser.parse_args()

    # Set up directories
    output_dir = Path(args.output_dir)
    label_directory = Path(args.input_file)

    # Determine if method uses DM embeddings
    # Check for n_dm argument and force_pca flag
    n_dm = getattr(args, "n_dm", 0)
    force_pca = getattr(args, "force_pca", False)
    use_dm = (n_dm > 0) and not force_pca

    # Load dataset with appropriate embeddings
    logger.info(f"Loading dataset for {args.method}")
    adata = data_loader.load_dataset(
        args.file_path,
        label_directory,
        args.ds_type,
        args.pop,
        args.pop_enr,
        args.seed,
        args.batch_sd,
        args.layer_embedding,
        use_dm=use_dm,
    )

    # Run method
    logger.step(f"Running {args.method}")
    results = method_module.run(adata, args)

    # Save main results
    logger.info(f"Saving results to {output_dir}")
    data_loader.save_results(
        results["results"],
        output_dir,
        args.ds_type,
        args.pop,
        args.pop_enr,
        args.seed,
        args.batch_sd,
        args.package,
        "_package_performance",
    )

    # Save any additional outputs
    if "additional" in results:
        for name, df in results["additional"].items():
            df.to_csv(output_dir / f"{name}.csv")

    logger.success(f"{args.method} completed successfully")


if __name__ == "__main__":
    main()
