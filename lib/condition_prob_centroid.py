import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.preprocessing import StandardScaler
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib.helper_functions import scale


def generate_enr_prob(pop, pop_enr, w_logit, cell_type_dict=None):
    """
    Assign enrichment score to the selected random_cells

    """

    n_clusters = w_logit.shape[1]
    enr_scores = pd.Series(0.5, index=w_logit.columns)

    if not isinstance(pop_enr, list):
        pop_enr_new = [pop_enr]
    else:
        pop_enr_new = pop_enr

    if not isinstance(pop, list):
        pop_new = [pop]
    else:
        pop_new = pop

    if len(pop_enr_new) == len(pop_new):
        # If 'pop_enr' is a list and matches the length of 'pop', assign directly
        for cluster, enr in zip(pop_new, pop_enr_new):
            enr_scores[cluster] = enr
    else:
        # If 'pop_enr' is a single value or doesn't match 'pop' in length, repeat 'pop_enr' for each 'pop' and assign
        # This ensures 'pop_enr' is treated as a repeated value if it's not already a list matching 'pop' in length
        pop_enr = np.repeat(pop_enr, len(pop_new))
        enr_scores[pop] = pop_enr

    return enr_scores


def normalize_enr_prob(w_logit, enr_scores):
    """
    Normalize enrichment probabilities for each cell population.

    Scales weight matrix values to range [0.5, enrichment_score] for each population.
    This creates smooth probability gradients where:
    - Cells far from centroid → probability = 0.5 (neutral)
    - Cells close to centroid → probability = enrichment_score

    Parameters:
    - w_logit: Weight matrix (cells x populations) from distance calculations
    - enr_scores: Enrichment scores for each population (0.5 = neutral, 0.95 = highly enriched)

    Returns:
    - enr_prob: Normalized probabilities (DataFrame, cells x populations)
    """
    enr_prob = pd.DataFrame(index=w_logit.index, columns=w_logit.columns)

    for i, col in enumerate(w_logit.columns):
        min_val = 0.5
        max_val = enr_scores.iloc[i]
        enr_prob[col] = scale(w_logit[col], min_val, max_val)
    return enr_prob


def set_relevant_prob(enr_prob, pop_enr, pop, adata, pop_column, cell_type_dict=None):
    """
    set probability to each cells when certain cell type is selected
    """

    # Initialize `cond_probability` with a default of 0.5
    prob_matrix = enr_prob[pop]

    # Initialize `cond_probability` with a default of 0.5
    cond_probability = pd.Series(0.5, index=enr_prob.index)

    if not isinstance(pop_enr, list):
        pop_enr_new = [pop_enr]
    else:
        pop_enr_new = pop_enr

    if not isinstance(pop, list):
        pop_new = [pop]
    else:
        pop_new = pop

        # If prob_matrix is not reduced to a single column, calculate the row means as the condition probability
    if len(pop_new) > 1:
        cond_probability = prob_matrix.mean(axis=1)

        # Update `cond_probability` for cells belonging to any of the populations in `pop`
        for population in pop_new:
            # Identify cells belonging to the current population
            cells_in_pop = adata.obs_names[adata.obs[pop_column] == population]

            # Directly assign probabilities from `prob_matrix` for these cells
            cond_probability.loc[cells_in_pop] = prob_matrix.loc[
                cells_in_pop, population
            ]
    else:
        cond_probability = prob_matrix
    return cond_probability
