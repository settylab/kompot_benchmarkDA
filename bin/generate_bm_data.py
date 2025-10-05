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

import argparse


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

    adata = read_file.read_dataset(args.file_path, args.layer_embedding)
    output_dir = Path(args.output_dir)

    if args.ds_type != "cluster":
        logger.debug(f"Processing non-cluster dataset: {args.ds_type}")
        X_emb = read_file.get_embedding_value(adata, args.layer_embedding)

        logger.debug("Computing weight matrix from centroids")
        w_logit, conditions = get_weight_matrix.get_weight_matrix_centroid(
            adata, args.pop_column, args.seed, X_emb, args.n_conditions, args.m, args.a_logit
        )

        logger.debug("Generating enrichment scores")
        enr_scores = condition_prob_centroid.create_enrichment_scores(
            args.pop, args.pop_enr, w_logit.columns
        )

        logger.debug("Computing condition probabilities")
        cond_probability = condition_prob_centroid.set_relevant_prob(
            w_logit, enr_scores, args.pop, adata, args.pop_column
        )

        logger.debug("Saving weight matrix")
        w_logit.to_csv(output_dir / "weight_matrix.csv")

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

        adata = synth_labels.quantile_assign_label_old(
            adata, args.pop_column, args.pop_enr, args.pop
        )

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

        X_pca = pd.DataFrame(
            adata.obsm[f"{args.layer_embedding}_batch"], index=adata.obs_names
        )

        # Assuming `adata` is your AnnData object
        obs_df = adata.obs.copy()

        # Add the index as a column named "rowname"
        obs_df["rowname"] = obs_df.index

        cols = ["rowname"] + [col for col in obs_df.columns if col != "rowname"]
        obs_df = obs_df[cols]

        obs_df.to_csv(
            output_dir
            / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}.coldata.csv",
            index=False,
        )

        str_batch = str(args.batch_sd)
        int_batch = helper_functions.convert_number_str(str_batch)

        # Save X_pca_batch (for PCA-based methods like MELD, CNA, R methods)
        X_pca["rowname"] = obs_df.index
        cols = ["rowname"] + [col for col in X_pca.columns if col != "rowname"]
        X_pca = X_pca[cols]
        X_pca.to_csv(
            output_dir
            / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}_batchEffect{int_batch}.emb.csv",
            index=False,
            float_format="%.16f",
        )

        # Save DM_EigenVectors_batch separately (for DM-based methods like Mellon, Kompot)
        if args.n_dm > 0:
            X_dm = pd.DataFrame(
                adata.obsm["DM_EigenVectors_batch"], index=adata.obs_names
            )
            X_dm["rowname"] = obs_df.index
            cols = ["rowname"] + [col for col in X_dm.columns if col != "rowname"]
            X_dm = X_dm[cols]
            X_dm.to_csv(
                output_dir
                / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}_batchEffect{int_batch}.emb.dm.csv",
                index=False,
                float_format="%.16f",
            )

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

        adata = cluster_dataset_synth_labels.quantile_assign_label_old(
            adata, args.pop_column, args.pop_enr, args.pop
        )

        logger.debug("Adding batch effects to embeddings")
        adata = cluster_dataset_synth_labels.add_batch_effect_pca(
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

        X_pca = pd.DataFrame(
            adata.obsm[f"{args.layer_embedding}_batch"], index=adata.obs_names
        )

        obs_df = adata.obs.copy()

        # Add the index as a column named "rowname"
        obs_df["rowname"] = obs_df.index
        cols = ["rowname"] + [col for col in obs_df.columns if col != "rowname"]
        obs_df = obs_df[cols]
        obs_df.to_csv(
            output_dir
            / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}.coldata.csv",
            index=False,
        )
        str_batch = str(args.batch_sd)
        int_batch = helper_functions.convert_number_str(str_batch)

        # Save X_pca_batch (for PCA-based methods like MELD, CNA, R methods)
        X_pca["rowname"] = obs_df.index
        cols = ["rowname"] + [col for col in X_pca.columns if col != "rowname"]
        X_pca = X_pca[cols]
        X_pca.to_csv(
            output_dir
            / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}_batchEffect{int_batch}.emb.csv",
            index=False,
            float_format="%.16f",
        )

        # Save DM_EigenVectors_batch separately (for DM-based methods like Mellon, Kompot)
        if args.n_dm > 0:
            X_dm = pd.DataFrame(
                adata.obsm["DM_EigenVectors_batch"], index=adata.obs_names
            )
            X_dm["rowname"] = obs_df.index
            cols = ["rowname"] + [col for col in X_dm.columns if col != "rowname"]
            X_dm = X_dm[cols]
            X_dm.to_csv(
                output_dir
                / f"benchmark_{args.ds_type}_pop_{args.pop}_enr{args.pop_enr}_seed{args.seed}_batchEffect{int_batch}.emb.dm.csv",
                index=False,
                float_format="%.16f",
            )

        logger.success(
            f"Successfully generated labels for {args.ds_type}-{args.pop}-{args.pop_enr}-{args.seed}-{args.batch_sd}"
        )


if __name__ == "__main__":
    main()
