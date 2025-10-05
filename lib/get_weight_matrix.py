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

from lib.find_centroid import find_centroid
from lib.euclidean_distance import euclidean_distance
from lib.weight_calculation import calculate_weights_centroid
from lib.helper_functions import log_function
from lib.constants import FUZZY_CMEANS_M, SIGMOID_STEEPNESS


def get_weight_matrix_centroid(
    adata, pop_column, seed, X_emb, n_conditions, m=FUZZY_CMEANS_M, a_logit=SIGMOID_STEEPNESS
):
    """
    Calculate sigmoid-transformed fuzzy membership weights from cell-to-centroid distances.

    This function implements a multi-step transformation pipeline:
    1. Find population centroids in embedding space
    2. Calculate Euclidean distances from cells to centroids
    3. Compute fuzzy membership weights using fuzzy c-means
    4. Normalize weights (z-score standardization)
    5. Apply sigmoid transformation to create smooth gradients

    The resulting weights represent each cell's "soft" membership to each population,
    transformed to have smooth sigmoid-shaped gradients suitable for enrichment modeling.

    Parameters:
    -----------
    adata : AnnData
        AnnData object containing cell metadata
    pop_column : str
        Column name in adata.obs containing population/cluster labels
    seed : int
        Random seed for reproducibility
    X_emb : DataFrame
        Cell embeddings (e.g., PCA coordinates), shape (n_cells, n_components)
    n_conditions : int
        Number of experimental conditions to generate
    m : float, default=FUZZY_CMEANS_M
        Fuzziness parameter for fuzzy c-means (m > 1)
    a_logit : float, default=SIGMOID_STEEPNESS
        Steepness parameter for sigmoid transformation

    Returns:
    --------
    sigmoid_fuzzy_weights : DataFrame
        Sigmoid-transformed fuzzy membership weights, shape (n_cells, n_populations)
        Values in (0, 1) range with smooth gradients
    conditions : list
        List of condition names ["Condition1", "Condition2", ...]

    Notes:
    ------
    Transformation pipeline:
    - fuzzy_weights: Raw fuzzy c-means weights (rows sum to ~1.0)
    - normalized_fuzzy_weights: Z-score normalized (mean=0, std=1 per column)
    - sigmoid_fuzzy_weights: Sigmoid transformed (smooth gradients in (0,1))
    """
    np.random.seed(seed)
    conditions = [f"Condition{i}" for i in range(1, n_conditions + 1)]

    # Step 1: Find population centroids in embedding space
    cluster_membership = adata.obs[pop_column]
    centroid_emb, cluster_names = find_centroid(X_emb, cluster_membership)

    # Step 2: Calculate Euclidean distances from cells to centroids
    centroid_distance = euclidean_distance(X_emb, centroid_emb)

    # Step 3: Compute fuzzy membership weights using fuzzy c-means
    # fuzzy_weights[i,j] = membership of cell i to population j
    fuzzy_weights = calculate_weights_centroid(centroid_distance, m)

    # Convert to DataFrame with proper labels
    fuzzy_weights_df = pd.DataFrame(
        fuzzy_weights,
        columns=cluster_names,
        index=X_emb.index
    )

    # Step 4: Normalize weights using z-score standardization
    # This ensures each population column has mean=0, std=1
    scaler = StandardScaler()
    normalized_fuzzy_weights = scaler.fit_transform(fuzzy_weights)

    # Step 5: Apply sigmoid transformation for smooth gradients
    # sigmoid(x) = 1 / (1 + exp(-a*x))
    # This maps normalized weights to (0, 1) with smooth S-curves
    sigmoid_fuzzy_weights = pd.DataFrame(
        log_function(normalized_fuzzy_weights, a=a_logit),
        columns=fuzzy_weights_df.columns,
        index=fuzzy_weights_df.index
    )

    return sigmoid_fuzzy_weights, conditions
