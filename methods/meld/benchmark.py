"""
MELD benchmark interface.
Contains method-specific argument parsing and execution logic.
"""

import pandas as pd
from . import runMELD


def add_arguments(parser):
    """Add MELD-specific arguments to parser."""
    parser.add_argument("--beta", type=float, help="Beta value")
    parser.add_argument("--k_meld", type=int, help="K MELD value")
    parser.add_argument(
        "--n_dm",
        type=int,
        default=0,
        help="Number of diffusion component for MELD (0 = PCA mode)",
    )


def run(adata, args):
    """
    Run MELD method and return results.

    Args:
        adata: AnnData object loaded with appropriate embeddings
        args: Parsed command-line arguments

    Returns:
        dict: Dictionary with 'results' DataFrame and optional additional outputs
    """
    # Determine embedding mode: DM by default, PCA if n_dm=0
    use_dm = args.n_dm > 0

    # Run MELD with hardcoded benchmark columns
    # synth_samples and synth_labels are a matched pair from the labeling process
    sample_likelihoods_meld, samplem = runMELD.runMELD(
        adata,
        args.k_meld,
        "synth_samples",  # Hardcoded - required for benchmark
        "synth_labels",  # Hardcoded - required for benchmark
        args.layer_embedding,
        args.beta,
        use_dm=use_dm,
        dm_comp=args.n_dm,
    )

    # Prepare results
    df_meld = pd.DataFrame(
        sample_likelihoods_meld,
        columns=[
            f"col_{i}" for i in range(sample_likelihoods_meld.reshape(-1, 1).shape[1])
        ],
        index=adata.obs_names,
    )

    return {"results": df_meld}
