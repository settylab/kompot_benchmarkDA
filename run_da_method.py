#!/usr/bin/env python3
"""
Simple CLI for running individual DA methods on datasets.
Provides an easy interface for modular execution without complex script generation.
"""

import argparse
import sys
import subprocess
from pathlib import Path
import os

# Add project root to Python path
PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.append(str(PROJECT_ROOT))

from config.method_config import PYTHON_METHODS, R_METHODS, list_methods_with_descriptions
from config.dataset_config import get_n_dm_for_dataset

def print_colored(text, color):
    colors = {
        'red': '\033[0;31m',
        'green': '\033[0;32m',
        'yellow': '\033[1;33m',
        'blue': '\033[0;34m',
        'bold': '\033[1m',
        'reset': '\033[0m'
    }
    print(f"{colors.get(color, '')}{text}{colors['reset']}")

def get_available_datasets():
    """Get list of available datasets."""
    datasets = []

    # Check synthetic datasets
    synthetic_datasets = ['linear', 'branch', 'cluster']
    for dataset in synthetic_datasets:
        synthetic_path = PROJECT_ROOT / 'data' / 'synthetic' / dataset / f'{dataset}.h5ad'
        if synthetic_path.exists():
            datasets.append(dataset)

    # Check real datasets
    real_datasets = ['bcr-xl', 'covid19-pbmc', 'levine32', 'pancreas']
    for dataset in real_datasets:
        real_path = PROJECT_ROOT / 'data' / 'real' / dataset / f'{dataset}.h5ad'
        if real_path.exists():
            datasets.append(dataset)

    return datasets

def validate_inputs(method, dataset, embedding_type):
    """Validate method, dataset, and embedding inputs."""
    errors = []

    # Check method
    if method not in PYTHON_METHODS and method not in R_METHODS:
        errors.append(f"Unknown method '{method}'. Available: {', '.join(list(PYTHON_METHODS.keys()) + list(R_METHODS.keys()))}")

    # Check dataset availability
    available_datasets = get_available_datasets()
    if dataset not in available_datasets:
        errors.append(f"Dataset '{dataset}' not available. Available: {', '.join(available_datasets)}")
        if not available_datasets:
            errors.append("No datasets found. Run preprocessing first: ./main_improved.sh preprocess")

    # Check embedding type
    if embedding_type not in ['pca', 'dm']:
        errors.append(f"Invalid embedding type '{embedding_type}'. Available: pca, dm")

    return errors

def run_python_method(method, dataset, embedding_type, **kwargs):
    """Run a Python DA method."""

    # Validate inputs
    validation_errors = validate_inputs(method, dataset, embedding_type)
    if validation_errors:
        for error in validation_errors:
            print_colored(f"ERROR: {error}", 'red')
        return False

    if method not in PYTHON_METHODS:
        print_colored(f"ERROR: '{method}' is not a Python method", 'red')
        print_colored(f"Available Python methods: {', '.join(PYTHON_METHODS.keys())}", 'blue')
        return False

    # Determine file path and parameters using dataset configuration
    n_dm = get_n_dm_for_dataset(dataset) if embedding_type == 'dm' else 0

    if dataset in ['covid19-pbmc', 'pancreas', 'bcr-xl', 'levine32']:
        data_file = PROJECT_ROOT / 'data' / 'real' / dataset / f'{dataset}_{embedding_type.upper()}_{n_dm}.h5ad'
    else:  # synthetic
        data_file = PROJECT_ROOT / 'data' / 'synthetic' / dataset / f'{dataset}_{embedding_type.upper()}_{n_dm}.h5ad'

    if not data_file.exists():
        print_colored(f"ERROR: Dataset file not found: {data_file}", 'red')
        print_colored("Run preprocessing first:", 'yellow')
        print_colored(f"  ./cli.sh --datasets {dataset} --embeddings {embedding_type} preprocess", 'yellow')
        return False

    # Create output directory
    output_dir = PROJECT_ROOT / 'results' / 'individual' / f'{method}_{dataset}_{embedding_type}'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Run the method
    script_path = PROJECT_ROOT / 'python_method' / PYTHON_METHODS[method]['script']

    cmd = [
        'python', str(script_path),
        '--file_path', str(data_file),
        '--pop', kwargs.get('population', 'M1'),
        '--pop_enr', str(kwargs.get('enrichment', 2.0)),
        '--pop_column', kwargs.get('pop_column', 'celltype'),
        '--ds_type', dataset,
        '--batch_sd', str(kwargs.get('batch_sd', 1.0)),
        '--input_file', str(data_file.parent / 'synthetic_labels'),  # placeholder
        '--package', method,
        '--seed', str(kwargs.get('seed', 0)),
        '--layer_embedding', 'X_pca',
        '--output_dir', str(output_dir)
    ]

    # Add method-specific parameters - always add n_dm for methods that need it
    if method in ['mellon', 'meld']:
        cmd.extend(['--n_dm', str(n_dm)])

    print_colored(f"Running {method} on {dataset} with {embedding_type} embedding...", 'blue')
    print_colored(f"Command: {' '.join(cmd)}", 'yellow')

    try:
        result = subprocess.run(cmd, cwd=PROJECT_ROOT, capture_output=True, text=True)
        if result.returncode == 0:
            print_colored(f"SUCCESS: Results saved to {output_dir}", 'green')
            return True
        else:
            print_colored(f"ERROR: Method failed", 'red')
            print_colored(f"STDOUT: {result.stdout}", 'yellow')
            print_colored(f"STDERR: {result.stderr}", 'red')
            return False
    except Exception as e:
        print_colored(f"ERROR: Failed to run method: {e}", 'red')
        return False

def run_r_method(method, dataset, embedding_type, **kwargs):
    """Run an R DA method."""

    # Validate inputs
    validation_errors = validate_inputs(method, dataset, embedding_type)
    if validation_errors:
        for error in validation_errors:
            print_colored(f"ERROR: {error}", 'red')
        return False

    if method not in R_METHODS:
        print_colored(f"ERROR: '{method}' is not an R method", 'red')
        print_colored(f"Available R methods: {', '.join(R_METHODS.keys())}", 'blue')
        return False

    # Check if Rscript is available
    import shutil
    if not shutil.which('Rscript'):
        print_colored("ERROR: Rscript not found. R methods require R to be installed.", 'red')
        print_colored("Install R or use the benchmarkda environment:", 'yellow')
        print_colored("  ./cli.sh setup", 'yellow')
        return False

    # Determine file paths - same as Python methods
    n_dm = get_n_dm_for_dataset(dataset) if embedding_type == 'dm' else 0

    if dataset in ['covid19-pbmc', 'pancreas', 'bcr-xl', 'levine32']:
        data_file = PROJECT_ROOT / 'data' / 'real' / dataset / f'{dataset}_{embedding_type.upper()}_{n_dm}.h5ad'
    else:  # synthetic
        data_file = PROJECT_ROOT / 'data' / 'synthetic' / dataset / f'{dataset}_{embedding_type.upper()}_{n_dm}.h5ad'

    if not data_file.exists():
        print_colored(f"ERROR: Dataset file not found: {data_file}", 'red')
        print_colored("Run preprocessing first:", 'yellow')
        print_colored(f"  ./cli.sh --datasets {dataset} --embeddings {embedding_type} preprocess", 'yellow')
        return False

    # Create output directory
    output_dir = PROJECT_ROOT / 'results' / 'individual' / f'{method}_{dataset}_{embedding_type}'
    output_dir.mkdir(parents=True, exist_ok=True)

    print_colored(f"Running R method '{method}' on {dataset} with {embedding_type} embedding", 'blue')
    print_colored("Note: R methods are best run through the full pipeline:", 'yellow')
    print_colored(f"  ./main_improved.sh --datasets {dataset} --methods r --embeddings {embedding_type} benchmark", 'yellow')

    return True

def main():
    parser = argparse.ArgumentParser(
        description='Run individual DA methods on datasets',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run kompot on linear dataset with DM embedding
  python run_da_method.py kompot linear --embedding dm

  # Run mellon on branch dataset with PCA embedding
  python run_da_method.py mellon branch --embedding pca --population M2

  # List available methods and datasets
  python run_da_method.py --list
        """
    )

    parser.add_argument('method', nargs='?', help='DA method to run')
    parser.add_argument('dataset', nargs='?', help='Dataset to use')
    parser.add_argument('--embedding', choices=['pca', 'dm'], default='dm',
                        help='Embedding type (default: dm)')
    parser.add_argument('--population', default='M1',
                        help='Population to test (default: M1)')
    parser.add_argument('--enrichment', type=float, default=2.0,
                        help='Enrichment level (default: 2.0)')
    parser.add_argument('--batch_sd', type=float, default=1.0,
                        help='Batch effect strength (default: 1.0)')
    parser.add_argument('--seed', type=int, default=0,
                        help='Random seed (default: 0)')
    parser.add_argument('--list', action='store_true',
                        help='List available methods and datasets')

    args = parser.parse_args()

    if args.list:
        print_colored("AVAILABLE METHODS WITH DESCRIPTIONS:", 'bold')
        print(list_methods_with_descriptions())

        print_colored("\nAVAILABLE DATASETS:", 'bold')
        datasets = get_available_datasets()
        if datasets:
            print_colored(f"Available: {', '.join(datasets)}", 'green')
        else:
            print_colored("No datasets found. Run preprocessing first:", 'yellow')
            print_colored("  ./cli.sh preprocess", 'yellow')

        print_colored("\nAVAILABLE EMBEDDINGS:", 'bold')
        print_colored("dm (diffusion maps), pca (principal components)", 'green')

        print_colored("\nEXAMPLES:", 'bold')
        if datasets:
            example_dataset = datasets[0]
            print_colored(f"  python run_da_method.py kompot {example_dataset} --embedding dm", 'blue')
            print_colored(f"  python run_da_method.py mellon {example_dataset} --embedding pca --population M2", 'blue')
        print_colored("  python run_da_method.py milo linear --embedding dm", 'blue')
        return

    if not args.method or not args.dataset:
        parser.print_help()
        return

    # Check if dataset exists
    available_datasets = get_available_datasets()
    if args.dataset not in available_datasets:
        print_colored(f"ERROR: Dataset '{args.dataset}' not found", 'red')
        print_colored(f"Available: {', '.join(available_datasets)}", 'blue')
        return

    # Run the method
    kwargs = {
        'population': args.population,
        'enrichment': args.enrichment,
        'batch_sd': args.batch_sd,
        'seed': args.seed,
        'pop_column': 'celltype'
    }

    if args.method in PYTHON_METHODS:
        success = run_python_method(args.method, args.dataset, args.embedding, **kwargs)
    elif args.method in R_METHODS:
        success = run_r_method(args.method, args.dataset, args.embedding, **kwargs)
    else:
        print_colored(f"ERROR: Unknown method '{args.method}'", 'red')
        print_colored("Use --list to see available methods", 'yellow')
        success = False

    if not success:
        sys.exit(1)

if __name__ == '__main__':
    main()