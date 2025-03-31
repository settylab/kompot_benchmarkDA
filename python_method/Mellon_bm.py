"""
Mellon method benchmark implementation.
Uses a standardized interface for running the Mellon method.
"""

import argparse
import numpy as np
import pandas as pd
from pathlib import Path

# Local imports
import runMellon
import data_loader
import helper_functions

def main():
    """Main function to run Mellon on a dataset."""
    # Create the parser
    parser = argparse.ArgumentParser(description='Run Mellon on a dataset for DA benchmark.')

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
    parser.add_argument('--n_dm', type=int, help='Number of diffusion component for Mellon value')
    parser.add_argument('--mellon_d_method', type=str, help='The d_method for Mellon, can be "embedding" or "fractal"')
    parser.add_argument('--norm_density', type=str, help='When using d_method as fractal, use normalized density or unnormalized')
    parser.add_argument('--hyperparameter', type=str, help='Whether use computed hyperparameter for Mellon')
    parser.add_argument('--corrected', type=str, help='Whether correct the log fold change mean of Mellon')
    parser.add_argument('--ls_factor', type=float, help='LS factor')
    parser.add_argument('--ls_mode', type=str, help='LS mode, PCA or DM')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')

    # Parse arguments
    args = parser.parse_args()

    # Set up directories
    output_dir = Path(args.output_dir)
    input_file = Path(args.input_file)

    # Process each iteration
    for i in range(1):
        iteration_directory = input_file / f'iteration_{i}'
        output_dir_i = output_dir / f'iteration_{i}'
        # Load dataset
        adata = data_loader.load_dataset(
            args.file_path, 
            iteration_directory, 
            args.ds_type, 
            args.pop, 
            args.pop_enr, 
            args.seed, 
            args.batch_sd, 
            args.layer_embedding
        )
        
        # Run Mellon with appropriate parameters
        if args.hyperparameter == "No":
            log_fold_change_mean, zscores = runMellon.runMELLON(
                adata, args.mellon_d_method, args.norm_density, "synth_labels", args.ls_mode, args.layer_embedding, args.n_dm, args.ls_factor
            )
        elif args.hyperparameter == "Yes":
            log_fold_change_mean, zscores = runMellon.runMELLON_synchronized(
                adata, args.mellon_d_method, args.norm_density, args.corrected, "synth_labels", args.n_dm, 
                args.ls_factor, args.ls_mode, args.layer_embedding
            )
        
        # Prepare results
        df_mellon_lfc = pd.DataFrame(
            log_fold_change_mean, 
            columns=[f"col_{i}" for i in range(log_fold_change_mean.reshape(-1,1).shape[1])],
            index=adata.obs_names
        )
        
        df_mellon_zscore = pd.DataFrame(
            zscores, 
            columns=[f"col_{i}" for i in range(zscores.reshape(-1,1).shape[1])],
            index=adata.obs_names
        )
        
        # Save results
        data_loader.save_results(
            df_mellon_lfc, 
            output_dir_i, 
            args.ds_type, 
            args.pop, 
            args.pop_enr, 
            args.seed, 
            args.batch_sd, 
            args.package, 
            "_package_performance"
        )
        
        #df_mellon_zscore.to_csv(output_dir_i / "mellon_zscore.csv")

if __name__ == "__main__":
    main()