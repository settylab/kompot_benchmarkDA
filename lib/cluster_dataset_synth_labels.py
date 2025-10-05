import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import palantir
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib.constants import SEED_BATCH_EFFECT, QUANTILE_LABEL_THRESHOLD_PCT


def add_synth_label_cluster_labels(
    adata, pop, seed, pop_enr, pop_col, n_conditions, cap_enr=None
):
    """
    Generate synthetic labels for cluster datasets.

    For cluster datasets, enrichment probabilities are assigned directly based on
    cluster membership (no distance calculations needed).

    Parameters:
    - adata: AnnData object with cluster labels in .obs[pop_col]
    - pop: Target population/cluster name(s) (string or list)
    - seed: Random seed for reproducibility
    - pop_enr: Enrichment probability for target population(s) (float or list)
    - pop_col: Column name containing cluster labels
    - n_conditions: Number of conditions (typically 2)
    - cap_enr: Optional maximum enrichment probability

    Returns:
    - adata: Modified AnnData with synth_labels, Condition1_prob, Condition2_prob in .obs
    """
    conditions = [f"Condition{i}" for i in range(1, n_conditions + 1)]
    enr_scores = pd.Series(0.5, index=adata.obs_names)

    if not isinstance(pop_enr, list):
        pop_enr_new = [pop_enr]
    else:
        pop_enr_new = pop_enr

    if not isinstance(pop, list):
        pop_new = [pop]
    else:
        pop_new = pop

    for cluster, enr in zip(pop_new, pop_enr_new):
        adata_temp = adata[adata.obs[pop_col] == cluster]
        cells_temp = list(adata_temp.obs_names)
        enr_scores[cells_temp] = enr

    cond_probability = enr_scores.copy()

    if cap_enr is not None:
        cond_probability = np.where(
            cond_probability > cap_enr, cap_enr, cond_probability
        )

    condition1_prob = cond_probability
    condition2_prob = 1 - condition1_prob

    cond_probability_df = pd.DataFrame(
        {conditions[0]: condition1_prob, conditions[1]: condition2_prob}
    )

    np.random.seed(seed)
    synth_labels = []

    for i in range(len(cond_probability_df)):
        probs = cond_probability_df.iloc[i].values
        label = np.random.choice(conditions, p=probs)
        synth_labels.append(label)

    adata.obs["synth_labels"] = synth_labels
    adata.obs["Condition1_prob"] = list(cond_probability_df.iloc[:, 0])
    adata.obs["Condition2_prob"] = list(cond_probability_df.iloc[:, 1])

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
    Assign true_labels based on Condition2 probability with fixed threshold.

    Uses a simple 10% offset from the enrichment level to create thresholds.
    This version is used for cluster datasets.

    Parameters:
    ----------
    adata : AnnData
        AnnData object with Condition2_prob in .obs
    pop_col : str
        Column name for population/cluster labels
    pop_enr : float
        Enrichment level for the target population (0 < pop_enr < 1)
    pop : str
        Target population name

    Returns:
    --------
    adata : AnnData
        Modified AnnData with true_labels in .obs ("NegLFC", "NotDA", or "PosLFC")
    """
    pop_tbl = pd.DataFrame(adata.obs[pop_col].value_counts())

    # Determine the quantile based on the condition
    if pop_enr < 0.5:
        da_lower = pop_enr + (pop_enr / 100) * QUANTILE_LABEL_THRESHOLD_PCT
        da_upper = 1 - da_lower
    else:
        da_upper = pop_enr - (pop_enr / 100) * QUANTILE_LABEL_THRESHOLD_PCT
        da_lower = 1 - da_upper
    print(da_lower, da_upper)
    assert da_upper > da_lower, "da_upper must be greater than da_lower"

    # Assign the synthetic labels based on 'Condition2_prob'
    true_label = np.where(
        np.array(adata.obs["Condition2_prob"]) < da_lower,
        "NegLFC",
        np.where(np.array(adata.obs["Condition2_prob"]) > da_upper, "PosLFC", "NotDA"),
    )

    adata.obs["true_labels"] = true_label
    # Replace spaces with underscores in pop if necessary
    pop = pop.replace(" ", "_")

    return adata


def add_batch_effect_pca(
    adata, layer_embedding, batch_col="synth_batches", norm_sd=0.5, seed=SEED_BATCH_EFFECT
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
