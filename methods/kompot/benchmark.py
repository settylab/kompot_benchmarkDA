"""
Kompot benchmark interface.
Contains method-specific argument parsing and execution logic.
"""

import pandas as pd
from pathlib import Path
from . import runKompot


def add_arguments(parser):
    """Add Kompot-specific arguments to parser."""
    parser.add_argument('--n_dm', type=int, default=10, help='Number of diffusion component for Kompot')
    parser.add_argument('--ls_factor', type=float, default=10.0, help='Length scale factor for Kompot')
    parser.add_argument('--n_landmarks', type=int, default=None, help='Number of landmarks for Kompot')
    parser.add_argument('--log_fold_change_threshold', type=float, default=1.0, help='Log fold change threshold')
    parser.add_argument('--pvalue_threshold', type=float, default=0.05, help='P-value threshold')
    parser.add_argument('--force_pca', action='store_true', help='Force PCA mode even if n_dm > 0')


def run(adata, args):
    """
    Run Kompot method and return results.

    Args:
        adata: AnnData object loaded with appropriate embeddings
        args: Parsed command-line arguments

    Returns:
        dict: Dictionary with 'results' DataFrame and optional additional outputs
    """
    # Determine embedding mode: DM by default, PCA if n_dm=0 or force_pca=True
    n_dm = args.n_dm if args.n_dm is not None else 10
    use_dm = (n_dm > 0) and not getattr(args, 'force_pca', False)

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

    return {
        'results': df_kompot_lfc,
        'additional': {'kompot_zscore': df_kompot_zscore}
    }
