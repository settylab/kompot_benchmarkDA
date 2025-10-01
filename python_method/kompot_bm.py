"""
Kompot method benchmark implementation.
Uses a standardized interface for running the Kompot method.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# Local imports
import runKompot
import data_loader
import helper_functions

def main():
    """Main function to run Kompot on a dataset."""
    # Create the parser
    parser = argparse.ArgumentParser(description='Run Kompot on a dataset for DA benchmark.')

    # Add arguments
    parser.add_argument('--file_path', type=str, help='Path to the file')
    parser.add_argument('--pop', type=str, help='Population')
    parser.add_argument('--pop_enr', type=float, help='Population enrichment')
    parser.add_argument('--pop_column', type=str, help='Population column')
    parser.add_argument('--ds_type', type=str, help='Type of dataset')
    parser.add_argument('--batch_sd', type=float, help='Batch standard deviation')
    parser.add_argument('--input_file', type=str, help='Input file path')
    parser.add_argument('--package', type=str, help='Which package is using')
    parser.add_argument('--seed', type=int, help='Seed for random number generation')
    parser.add_argument('--layer_embedding', type=str, help='Layer embedding, X_pca or DM_EigenVectors')
    parser.add_argument('--n_dm', type=int, default=10, help='Number of diffusion component for Kompot')
    parser.add_argument('--ls_factor', type=float, default=10.0, help='Length scale factor for Kompot')
    parser.add_argument('--n_landmarks', type=int, default=None, help='Number of landmarks for Kompot')
    parser.add_argument('--log_fold_change_threshold', type=float, default=1.0, help='Log fold change threshold')
    parser.add_argument('--pvalue_threshold', type=float, default=0.05, help='P-value threshold')
    parser.add_argument('--force_pca', action='store_true', help='Force PCA mode even if n_dm > 0')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')

    # Parse arguments
    args = parser.parse_args()

    # Set up directories
    output_dir = Path(args.output_dir)
    label_directory = Path(args.input_file)

    # Determine embedding mode: DM by default, PCA if n_dm=0 or force_pca=True
    n_dm = args.n_dm if args.n_dm is not None else 10
    use_dm = (n_dm > 0) and not getattr(args, 'force_pca', False)

    # Load dataset with appropriate embeddings
    adata = data_loader.load_dataset(
        args.file_path,
        label_directory,
        args.ds_type,
        args.pop,
        args.pop_enr,
        args.seed,
        args.batch_sd,
        args.layer_embedding,
        use_dm=use_dm
    )

    # Run Kompot with standardized embedding handling
    log_fold_change_mean, zscores = runKompot.runKOMPOT_with_params(
        adata=adata,
        label_col="synth_labels",
        use_dm=use_dm,
        dm_comp=n_dm,
        ls_factor=args.ls_factor,
        n_landmarks=args.n_landmarks,
        log_fold_change_threshold=args.log_fold_change_threshold,
        pvalue_threshold=args.pvalue_threshold,
        random_state=args.seed
    )

    # Prepare results
    df_kompot_lfc = pd.DataFrame(
        log_fold_change_mean.reshape(-1, 1),
        columns=[f"col_{i}" for i in range(log_fold_change_mean.reshape(-1,1).shape[1])],
        index=adata.obs_names
    )

    df_kompot_zscore = pd.DataFrame(
        zscores.reshape(-1, 1),
        columns=[f"col_{i}" for i in range(zscores.reshape(-1,1).shape[1])],
        index=adata.obs_names
    )

    # Save results
    data_loader.save_results(
        df_kompot_lfc,
        output_dir,
        args.ds_type,
        args.pop,
        args.pop_enr,
        args.seed,
        args.batch_sd,
        args.package,
        "_package_performance"
    )

    df_kompot_zscore.to_csv(output_dir / "kompot_zscore.csv")

if __name__ == "__main__":
    main()