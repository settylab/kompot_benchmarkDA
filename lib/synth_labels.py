import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import os
import sys


def cap_probabilities(adata, cond_probability, conditions, cap_enr=None):
    """
    Generate condition probabilities for each cell.

    Parameters:
    - adata: AnnData object
    - cond_probability: Array of condition 1 probabilities for each cell
    - conditions: List of condition names (e.g., ["Condition1", "Condition2"])
    - cap_enr: Optional maximum enrichment probability (capped if exceeded)

    Returns:
    - conditions: List of condition names
    - cond_probability_df: DataFrame with probabilities for both conditions
    """
    if cap_enr is not None:
        cond_probability = np.where(
            cond_probability > cap_enr, cap_enr, cond_probability
        )

    condition1_prob = cond_probability
    condition2_prob = 1 - condition1_prob

    cond_probability_df = pd.DataFrame(
        {conditions[0]: condition1_prob, conditions[1]: condition2_prob}
    )
    return conditions, cond_probability_df


def label_condition_and_rep_labels(adata, cond_probability, seed):
    """
    Assign synth_labels to each cell based on condition probabilities.

    Uses weighted random sampling where each cell is assigned to a condition
    based on its condition probabilities.

    Parameters:
    - adata: AnnData object
    - cond_probability: DataFrame with condition probabilities (columns = conditions)
    - seed: Random seed for reproducibility

    Returns:
    - adata: Modified AnnData with synth_labels, Condition1_prob, Condition2_prob in .obs
    """
    np.random.seed(seed)
    conditions = cond_probability.columns.tolist()
    synth_labels = []

    for i in range(len(cond_probability)):
        probs = cond_probability.iloc[i].values
        label = np.random.choice(conditions, p=probs)
        synth_labels.append(label)

    adata.obs["synth_labels"] = synth_labels
    adata.obs["Condition1_prob"] = list(cond_probability.iloc[:, 0])
    adata.obs["Condition2_prob"] = list(cond_probability.iloc[:, 1])

    return adata


def label_condition_and_rep_other(adata, n_replicates, n_batches, seed):
    """
    Generate sample and batch labels from existing synth_labels.

    Creates:
    - synth_samples: Combination of condition label and replicate (e.g., "Condition1_R1")
    - synth_batches: Batch assignments for each sample

    Parameters:
    - adata: AnnData object with synth_labels already assigned
    - n_replicates: Number of replicates per condition
    - n_batches: Number of batches
    - seed: Random seed for batch shuffling

    Returns:
    - adata: Modified AnnData with synth_samples and synth_batches in .obs
    """
    np.random.seed(seed)
    synth_labels = adata.obs["synth_labels"]

    replicates = [f"R{i}" for i in range(1, n_replicates + 1)]

    batch_labels = [
        f"B{i}" for i in range(1, n_batches + 1) for _ in range(n_replicates)
    ]
    np.random.shuffle(batch_labels)
    batches = batch_labels

    synth_samples = [
        f"{label}_{rep}"
        for label, rep in zip(synth_labels, replicates * len(synth_labels))
    ]

    if n_batches > 1:
        batch_labels_dict = {
            sample: batches[i] for i, sample in enumerate(sorted(set(synth_samples)))
        }
    else:
        batch_labels_dict = {sample: "B1" for sample in set(synth_samples)}
    synth_batches = [batch_labels_dict[sample] for sample in synth_samples]

    adata.obs["synth_samples"] = synth_samples
    adata.obs["synth_batches"] = synth_batches

    return adata


def quantile_assign_label(adata, pop_col, pop_enr, pop):
    """
    Assign true_labels for evaluation based on condition probabilities.

    Uses quantile-based thresholds to classify cells as:
    - "PosLFC": Enriched in condition 1
    - "NegLFC": Enriched in condition 2
    - "NotDA": Not differentially abundant

    Parameters:
    - adata: AnnData object with Condition1_prob
    - pop_col: Column name for population/cluster labels
    - pop_enr: Enrichment level for the target population
    - pop: Target population name

    Returns:
    - adata: Modified AnnData with true_labels in .obs
    """
    pop_tbl = adata.obs[pop_col].value_counts()

    total_cells = pop_tbl.sum()
    n_enr_cells = pop_tbl[pop]
    neutral_prop = (total_cells - n_enr_cells) / (total_cells * 2)
    print(neutral_prop)
    # Determine the quantile based on the condition
    if pop_enr < 0.5:
        da_lower = np.quantile(
            adata.obs["Condition1_prob"], pop_tbl[pop_tbl.index == pop] / total_cells
        )[0]
        da_upper = (2 * neutral_prop) - da_lower
    else:
        da_upper = np.quantile(
            adata.obs["Condition1_prob"],
            1 - pop_tbl[pop_tbl.index == pop] / total_cells,
        )[0]
        print(da_upper)

        da_lower = (2 * neutral_prop) - da_upper
        print(da_lower)

    assert da_upper > da_lower, "da_upper must be greater than da_lower"

    bins = [-float("inf"), da_lower, da_upper, float("inf")]
    labels = ["PosLFC", "NotDA", "NegLFC"]
    adata.obs["true_labels"] = pd.cut(adata.obs["Condition1_prob"], bins, labels=labels)

    return adata


def quantile_assign_label_old(adata, pop_col, pop_enr, pop):
    """
    DEPRECATED: Old version of quantile_assign_label.
    Use quantile_assign_label() instead (inverted probability logic).
    """
    pop_tbl = pd.DataFrame(adata.obs[pop_col].value_counts())

    # Determine the quantile based on the condition
    if pop_enr < 0.5:
        da_lower = np.quantile(
            adata.obs["Condition1_prob"],
            pop_tbl[pop_tbl.index == pop] / sum(pop_tbl.iloc[:, 0]),
        )
        da_upper = 1 - da_lower
    else:
        da_upper = np.quantile(
            adata.obs["Condition1_prob"],
            1 - pop_tbl[pop_tbl.index == pop] / sum(pop_tbl.iloc[:, 0]),
        )
        da_lower = 1 - da_upper

    assert da_upper > da_lower, "da_upper must be greater than da_lower"

    # Assign the synthetic labels based on 'Condition2_prob'
    true_label = np.where(
        np.array(adata.obs["Condition2_prob"]) < da_lower,
        "NegLFC",
        np.where(np.array(adata.obs["Condition2_prob"]) > da_upper, "PosLFC", "NotDA"),
    )

    adata.obs["true_labels"] = true_label[0]
    # Replace spaces with underscores in pop if necessary
    pop = pop.replace(" ", "_")

    return adata


def add_batch_effect_pca(
    adata, layer_embedding, batch_col="synth_batches", norm_sd=0.5, seed=43
):
    """
    Add batch effects to embeddings by adding batch-specific noise.

    Each batch gets a unique random offset scaled by the data variance.
    This simulates technical batch effects while preserving biological signal.

    Batch effects are scaled by sqrt(sum of component variances) to make
    norm_sd comparable across datasets with different scales.

    Parameters:
    - adata: AnnData object with embeddings in .obsm[layer_embedding]
    - layer_embedding: Name of embedding ('X_pca' or 'DM_EigenVectors')
    - batch_col: Column in .obs with batch labels (default: 'synth_batches')
    - norm_sd: Batch effect strength as fraction of data variance
               (0 = no effect, 1 = effect size equals data std, default: 0.5)
    - seed: Random seed for reproducibility

    Returns:
    - adata: Modified AnnData with batch-affected embeddings in .obsm['{layer_embedding}_batch']
    """
    np.random.seed(seed)
    # Extract embeddings
    X_emb = adata.obsm[layer_embedding]

    X_emb_df = pd.DataFrame(X_emb, index=adata.obs_names)

    # Initialize batch-corrected version with the original
    X_emb_batch = X_emb_df.copy()

    # Calculate scaling factor: sqrt of sum of variances across all components
    # This makes batch effect strength comparable across different datasets
    component_vars = np.var(X_emb, axis=0)
    variance_scale = np.sqrt(np.sum(component_vars))

    # Split cells by batch
    cell_batches = adata.obs[batch_col]

    # Add batch effect for each batch
    # norm_sd=0 means no effect, norm_sd=1 means effect ~ data std
    for batch in cell_batches.unique():
        # Generate batch effect scaled by data variance
        batch_effect = np.random.normal(
            loc=0, scale=norm_sd * variance_scale, size=X_emb.shape[1]
        )
        batch_indices = list(cell_batches[cell_batches == batch].index)
        X_emb_batch.loc[batch_indices] += batch_effect

    # Store the modified embeddings with batch effects
    adata.obsm[f"{layer_embedding}_batch"] = X_emb_batch.values

    return adata


### def add_batch_effect_DM(adata, layer_embedding,batch_col="synth_batches", norm_sd=0.5,seed = 43):
