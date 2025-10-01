"""
CNA method benchmark implementation.
Uses a standardized interface for running the CNA method.
"""

import argparse
import numpy as np
import pandas as pd
import anndata
from pathlib import Path

# Local imports
import runCNA
import data_loader
import helper_functions

def main():
    """Main function to run CNA on a dataset."""
    # Create the parser
    parser = argparse.ArgumentParser(description='Run CNA on a dataset for DA benchmark.')

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
    parser.add_argument('--k_cna', type=int, help='K CNA value')
    parser.add_argument('--n_dm', type=int, default=0, help='Number of diffusion component for CNA (0 = PCA mode)')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')

    # Parse arguments
    args = parser.parse_args()

    # Set up directories
    output_dir = Path(args.output_dir)
    label_directory = Path(args.input_file)

    # Load dataset (CNA is PCA-based, always use PCA embeddings)
    adata = data_loader.load_dataset(
        args.file_path,
        label_directory,
        args.ds_type,
        args.pop,
        args.pop_enr,
        args.seed,
        args.batch_sd,
        args.layer_embedding,
        use_dm=False
    )

    # Determine embedding mode: DM by default, PCA if n_dm=0
    use_dm = (args.n_dm > 0)

    # Run CNA with standardized embedding handling
    cna_res, md = runCNA.runCNA_func(
        adata,
        args.k_cna,
        "synth_samples",
        "synth_labels",
        {'Condition1': 0, "Condition2": 1},
        "synth_batches",
        layer_embedding=args.layer_embedding,
        use_dm=use_dm,
        dm_comp=args.n_dm
    )

    # Prepare results
    da_cell = np.repeat(0., len(md))
    da_cell[cna_res.kept] = cna_res.ncorrs

    df_cna_da = pd.DataFrame(
        da_cell,
        columns=[f"col_{i}" for i in range(da_cell.reshape(-1,1).shape[1])],
        index=adata.obs_names
    )

    # Save results
    data_loader.save_results(
        df_cna_da,
        output_dir,
        args.ds_type,
        args.pop,
        args.pop_enr,
        args.seed,
        args.batch_sd,
        args.package,
        "_package_performance"
    )

if __name__ == "__main__":
    main()