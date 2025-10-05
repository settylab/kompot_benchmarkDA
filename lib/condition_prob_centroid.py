"""
Condition probability calculation for synthetic differential abundance labeling.

This module handles the generation of enrichment probabilities for synthetic
benchmark data. It computes smooth transcriptional gradients based on distance
to population centroids, ensuring biologically realistic enrichment patterns.

Key functions:
- create_enrichment_scores: Map populations to target enrichment levels
- set_relevant_prob: Compute final probabilities with smooth gradients from fuzzy membership scores
"""
import numpy as np
import pandas as pd
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib.helper_functions import scale


def create_enrichment_scores(pop, pop_enr, populations):
    """
    Create enrichment score mapping for each population.

    Assigns target enrichment values to specified populations,
    with all other populations defaulting to 0.5 (neutral).

    Parameters:
    - pop: Population name(s) to enrich (string or list)
    - pop_enr: Enrichment level(s) for target populations (float or list)
    - populations: All population names (list or Index)

    Returns:
    - enr_scores: Series mapping population names to enrichment scores
                  (target populations → specified enrichment, others → 0.5)

    Example:
        create_enrichment_scores('M1', 0.95, ['M1', 'M2', 'M3'])
        → Series({'M1': 0.95, 'M2': 0.5, 'M3': 0.5})
    """
    # Initialize all populations to neutral (0.5)
    enr_scores = pd.Series(0.5, index=populations)

    # Ensure inputs are lists
    pop_list = [pop] if not isinstance(pop, list) else pop
    enr_list = [pop_enr] if not isinstance(pop_enr, list) else pop_enr

    # Assign enrichment scores to target populations
    if len(enr_list) == len(pop_list):
        # One enrichment value per population
        for population, enrichment in zip(pop_list, enr_list):
            enr_scores[population] = enrichment
    else:
        # Single enrichment value for all target populations
        for population in pop_list:
            enr_scores[population] = enr_list[0]

    return enr_scores


def set_relevant_prob(sigmoid_fuzzy_weights, enr_scores, pop, adata, pop_column):
    """
    Compute final enrichment probability for each cell based on proximity to enriched centroids.

    Uses a consistent algorithm for all cells to ensure smooth transcriptional gradients
    without artificial discontinuities at cluster label boundaries. This function combines
    the scaling and rescaling logic that was previously split between normalize_enr_prob
    and set_relevant_prob.

    Algorithm:
    1. Scale each enriched population column to [0, target_enrichment_score]
    2. Take max probability across enriched populations (minimum distance principle)
    3. Find minimum probability among cells labeled as enriched populations
    4. Rescale [min_prob, max_prob] → [0.5, max_prob] to ensure labeled cells ≥ 0.5
    5. Floor all probabilities at 0.5 (neutral baseline for distant cells)

    This logic applies consistently to both single and multiple enriched populations.

    Parameters:
    ----------
    sigmoid_fuzzy_weights : DataFrame
        Sigmoid-transformed fuzzy membership weights (cells x populations)
        from get_weight_matrix_centroid(). Values in (0, 1) range.
    enr_scores : Series
        Mapping of population names to target enrichment scores (0.5-1.0)
    pop : str or list
        Population name(s) to enrich
    adata : AnnData
        AnnData object with cluster labels in .obs[pop_column]
    pop_column : str
        Column name containing cluster labels

    Returns:
    --------
    cond_probability : Series
        Final enrichment probability for each cell

    Example:
        For single population M1 with enrichment 0.95:
        - Cells close to M1 centroid → ~0.95
        - Cells far from M1 centroid but labeled M1 → ~0.5
        - Other cells → 0.5
    """
    # Ensure pop is a list
    if not isinstance(pop, list):
        pop = [pop]

    # Select columns for enriched populations
    prob_matrix = sigmoid_fuzzy_weights[pop].copy()

    # Scale each column to [0, target_enrichment_score]
    for col in prob_matrix.columns:
        prob_matrix[col] = scale(prob_matrix[col], 0, enr_scores[col])

    # Take max probability (closest enriched centroid)
    cond_probability = prob_matrix.max(axis=1) if len(pop) > 1 else prob_matrix.iloc[:, 0]

    # Find minimum probability among cells labeled as enriched populations
    enriched_mask = adata.obs[pop_column].isin(pop)

    if enriched_mask.any():
        min_prob = cond_probability[enriched_mask].min()
        max_prob = cond_probability.max()

        # Rescale so minimum in enriched populations becomes 0.5
        # Linear rescaling: [min_prob, max_prob] → [0.5, max_prob]
        if max_prob > min_prob:
            cond_probability = 0.5 + (cond_probability - min_prob) / (max_prob - min_prob) * (max_prob - 0.5)

    # Floor at 0.5 (neutral baseline for cells far from enriched centroids)
    cond_probability = np.maximum(cond_probability, 0.5)

    return cond_probability
