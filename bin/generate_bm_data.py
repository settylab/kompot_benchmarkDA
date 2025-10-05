import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib import condition_prob_centroid
from lib import calculate_diffusion_map
from lib import get_weight_matrix
from lib import helper_functions
from lib import read_file
from lib.logger import get_logger
from lib.constants import (
    validate_enrichment,
    validate_m_parameter,
    validate_batch_sd,
    FUZZY_CMEANS_M,
    SIGMOID_STEEPNESS
)

import argparse


def save_benchmark_outputs(adata, args, output_dir, logger):
    """
    Save benchmark outputs: coldata and embeddings.

    Common function for both cluster and non-cluster datasets to save:
    - Cell metadata (coldata.csv)
    - PCA/batch-affected embeddings (.emb.csv)
    - DM embeddings if computed (.emb.dm.csv)

    Parameters:
    -----------
    adata : AnnData
        AnnData object with batch-affected embeddings in .obsm
    args : Namespace
        Parsed command-line arguments
    output_dir : Path
        Output directory for saving files
    logger : Logger
        Logger instance for progress messages

    Returns:
    --------
    None (files are saved to disk)
    """
    logger.debug("Saving coldata and embeddings")

    # Extract and save cell metadata
    obs_df = adata.obs.copy()
    obs_df["rowname"] = obs_df.index
    cols = ["rowname"] + [col for col in obs_df.columns if col != "rowname"]
    obs_df = obs_df[cols]

    coldata_path = output_dir / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}.coldata.csv"
    obs_df.to_csv(coldata_path, index=False)
    logger.debug(f"Saved coldata to {coldata_path.name}")

    # Convert batch_sd to int if it's a whole number (for filename)
    str_batch = str(args.batch_sd)
    int_batch = helper_functions.convert_number_str(str_batch)

    # Save PCA/batch-affected embeddings
    X_pca = pd.DataFrame(
        adata.obsm[f"{args.layer_embedding}_batch"], index=adata.obs_names
    )
    X_pca["rowname"] = obs_df.index
    cols = ["rowname"] + [col for col in X_pca.columns if col != "rowname"]
    X_pca = X_pca[cols]

    emb_path = output_dir / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}_batchEffect{int_batch}.emb.csv"
    X_pca.to_csv(emb_path, index=False, float_format="%.16f")
    logger.debug(f"Saved PCA embeddings to {emb_path.name}")

    # Save DM embeddings if computed
    if args.n_dm > 0:
        X_dm = pd.DataFrame(
            adata.obsm["DM_EigenVectors_batch"], index=adata.obs_names
        )
        X_dm["rowname"] = obs_df.index
        cols = ["rowname"] + [col for col in X_dm.columns if col != "rowname"]
        X_dm = X_dm[cols]

        dm_path = output_dir / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}_batchEffect{int_batch}.emb.dm.csv"
        X_dm.to_csv(dm_path, index=False, float_format="%.16f")
        logger.debug(f"Saved DM embeddings to {dm_path.name}")


def add_batch_effects_and_compute_dm(adata, args, logger):
    """
    Add batch effects to embeddings and optionally compute diffusion map.

    Common function for both cluster and non-cluster datasets to:
    1. Add batch effects to original embeddings
    2. Optionally compute DM on batch-affected embeddings

    Parameters:
    -----------
    adata : AnnData
        AnnData object with original embeddings
    args : Namespace
        Parsed command-line arguments
    logger : Logger
        Logger instance for progress messages

    Returns:
    --------
    adata : AnnData
        Modified AnnData with batch-affected embeddings (and DM if requested)
    """
    logger.debug("Adding batch effects to embeddings")
    adata = synth_labels.add_batch_effect_pca(
        adata,
        args.layer_embedding,
        batch_col="synth_batches",
        norm_sd=args.batch_sd,
        seed=args.seed,
    )

    # Compute DM on batch-simulated embeddings if needed
    if args.n_dm > 0:
        logger.info(
            f"Computing DM ({args.n_dm} components) on batch-simulated embeddings"
        )
        adata = calculate_diffusion_map.calculate_dm(
            adata, f"{args.layer_embedding}_batch", args.n_dm
        )
        logger.debug("DM computation complete")

    return adata


def main():
    # Create the parser
    parser = argparse.ArgumentParser(
        description="Generate benchmark data with synthetic condition labels and batch effects."
    )

    # Add arguments
    parser.add_argument("--file_path", type=str, help="Path to the file")

    parser.add_argument("--pop", type=str, help="Population")
    parser.add_argument("--pop_enr", type=float, help="Population enrichment")
    parser.add_argument("--pop_column", type=str, help="Population column")

    parser.add_argument("--ds_type", type=str, help="type of dataset")
    parser.add_argument("--batch_sd", type=float, help="batch standard deviation")

    parser.add_argument("--n_conditions", type=int, help="Number of conditions")
    parser.add_argument("--n_replicates", type=int, help="Number of replicates")
    parser.add_argument("--n_batches", type=int, help="Number of batches")
    parser.add_argument("--seed", type=int, help="Seed for random number generation")

    parser.add_argument("--m", type=float, help="M value")
    parser.add_argument("--a_logit", type=float, help="A_logit value")

    parser.add_argument(
        "--layer_embedding", type=str, help="Layer embedding, X_pca or DM_EigenVectors"
    )
    parser.add_argument(
        "--n_dm",
        type=int,
        default=0,
        help="Number of diffusion components to compute on batch-simulated embeddings",
    )
    parser.add_argument(
        "--output_dir", type=str, required=True, help="Output directory path"
    )

    # Parse arguments
    args = parser.parse_args()

    from lib import synth_labels
    from lib import cluster_dataset_synth_labels

    # Initialize logger
    logger = get_logger("generate_bm_data", verbose=False)

    logger.info(
        f"Dataset: {args.ds_type}, Population: {args.pop}, Enrichment: {args.pop_enr}, Seed: {args.seed}, Batch SD: {args.batch_sd}"
    )
    logger.debug(f"Output directory: {args.output_dir}")

    # ========================================================================
    # Input Validation
    # ========================================================================
    logger.debug("Validating input parameters")

    # Validate enrichment probability
    try:
        validate_enrichment(args.pop_enr)
    except ValueError as e:
        logger.error(f"Invalid enrichment value: {e}")
        sys.exit(1)

    # Validate batch SD
    try:
        validate_batch_sd(args.batch_sd)
    except ValueError as e:
        logger.error(f"Invalid batch SD: {e}")
        sys.exit(1)

    # Validate m parameter
    try:
        validate_m_parameter(args.m)
    except ValueError as e:
        logger.error(f"Invalid m parameter: {e}")
        sys.exit(1)

    # Validate file path exists
    file_path = Path(args.file_path)
    if not file_path.exists():
        logger.error(f"Data file not found: {file_path}")
        sys.exit(1)

    # Validate output directory can be created
    output_dir = Path(args.output_dir)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logger.error(f"Cannot create output directory {output_dir}: {e}")
        sys.exit(1)

    # Validate population is not empty
    if not args.pop or (isinstance(args.pop, str) and args.pop.strip() == ""):
        logger.error("Population name cannot be empty")
        sys.exit(1)

    # Validate number of conditions, replicates, batches
    if args.n_conditions < 2:
        logger.error(f"Number of conditions must be >= 2, got {args.n_conditions}")
        sys.exit(1)

    if args.n_replicates < 1:
        logger.error(f"Number of replicates must be >= 1, got {args.n_replicates}")
        sys.exit(1)

    if args.n_batches < 1:
        logger.error(f"Number of batches must be >= 1, got {args.n_batches}")
        sys.exit(1)

    if args.n_dm < 0:
        logger.error(f"Number of DM components must be >= 0, got {args.n_dm}")
        sys.exit(1)

    logger.debug("All input parameters validated successfully")

    adata = read_file.read_dataset(args.file_path, args.layer_embedding)
    output_dir = Path(args.output_dir)

    if args.ds_type != "cluster":
        logger.debug(f"Processing non-cluster dataset: {args.ds_type}")
        X_emb = read_file.get_embedding_value(adata, args.layer_embedding)

        logger.debug("Computing weight matrix from centroids")
        sigmoid_fuzzy_weights, conditions = get_weight_matrix.get_weight_matrix_centroid(
            adata, args.pop_column, args.seed, X_emb, args.n_conditions, args.m, args.a_logit
        )

        logger.debug("Generating enrichment scores")
        enr_scores = condition_prob_centroid.create_enrichment_scores(
            args.pop, args.pop_enr, sigmoid_fuzzy_weights.columns
        )

        logger.debug("Computing condition probabilities")
        cond_probability = condition_prob_centroid.set_relevant_prob(
            sigmoid_fuzzy_weights, enr_scores, args.pop, adata, args.pop_column
        )

        logger.step("Generating synthetic labels and batch-simulated embeddings")

        # Generate labels directly to output_dir (no iteration subdirectory)
        if not isinstance(cond_probability, pd.DataFrame):
            cond_probability_temp = pd.DataFrame(
                cond_probability, index=adata.obs_names
            )
        else:
            cond_probability_temp = cond_probability

        temp_cond = cond_probability_temp.iloc[:, 0]

        conditions, cond_probability_df = synth_labels.cap_probabilities(
            adata, temp_cond, conditions, cap_enr=None
        )

        adata = synth_labels.label_condition_and_rep_labels(
            adata, cond_probability_df, args.seed
        )

        adata = synth_labels.label_condition_and_rep_other(
            adata, args.n_replicates, args.n_batches, args.seed
        )

        adata = synth_labels.quantile_assign_label(
            adata, args.pop_column, args.pop_enr, args.pop
        )

        # Add batch effects and compute DM (common function)
        adata = add_batch_effects_and_compute_dm(adata, args, logger)

        # Save outputs (common function)
        save_benchmark_outputs(adata, args, output_dir, logger)

        logger.success(
            f"Successfully generated labels for {args.ds_type}-{args.pop}-{args.pop_enr}-{args.seed}-{args.batch_sd}"
        )

    elif args.ds_type == "cluster":
        logger.debug(f"Processing cluster dataset: {args.ds_type}")

        adata = cluster_dataset_synth_labels.add_synth_label_cluster_labels(
            adata, args.pop, args.seed, args.pop_enr, args.pop_column, args.n_conditions, cap_enr=None
        )
        adata = cluster_dataset_synth_labels.label_condition_and_rep_other(
            adata, args.n_replicates, args.n_batches, args.seed
        )

        adata = cluster_dataset_synth_labels.quantile_assign_label(
            adata, args.pop_column, args.pop_enr, args.pop
        )

        # Add batch effects and compute DM (common function)
        adata = add_batch_effects_and_compute_dm(adata, args, logger)

        # Save outputs (common function)
        save_benchmark_outputs(adata, args, output_dir, logger)

        logger.success(
            f"Successfully generated labels for {args.ds_type}-{args.pop}-{args.pop_enr}-{args.seed}-{args.batch_sd}"
        )


if __name__ == "__main__":
    main()
