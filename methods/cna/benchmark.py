"""
CNA benchmark interface.
Contains method-specific argument parsing and execution logic.
"""

import numpy as np
import pandas as pd
from . import runCNA


def add_arguments(parser):
    """Add CNA-specific arguments to parser."""
    parser.add_argument("--k_cna", type=int, help="K CNA value")
    parser.add_argument(
        "--n_dm",
        type=int,
        default=0,
        help="Number of diffusion component for CNA (0 = PCA mode)",
    )


def run(adata, args):
    """
    Run CNA method and return results.

    Args:
        adata: AnnData object loaded with appropriate embeddings
        args: Parsed command-line arguments

    Returns:
        dict: Dictionary with 'results' DataFrame and optional additional outputs
    """
    # Determine embedding mode: DM by default, PCA if n_dm=0
    use_dm = args.n_dm > 0

    # Run CNA with standardized embedding handling
    cna_res, md = runCNA.runCNA_func(
        adata,
        args.k_cna,
        "synth_samples",
        "synth_labels",
        {"Condition1": 0, "Condition2": 1},
        "synth_batches",
        layer_embedding=args.layer_embedding,
        use_dm=use_dm,
        dm_comp=args.n_dm,
    )

    # Prepare results
    da_cell = np.repeat(0.0, len(md))
    da_cell[cna_res.kept] = cna_res.ncorrs

    df_cna_da = pd.DataFrame(
        da_cell,
        columns=[f"col_{i}" for i in range(da_cell.reshape(-1, 1).shape[1])],
        index=adata.obs_names,
    )

    return {"results": df_cna_da}
