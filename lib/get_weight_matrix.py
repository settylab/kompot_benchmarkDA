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
    Calculated final weight matrix when the distance is calculated from cells to centroids
    1. adata: original adata to work on
    2. pop: which population
    3. pop_column: the columns for population in adata.obs
    4. seed: random seed
    5. X_emb: cells on embedding space
    6. mode: can be "euclidean" or "KNN", which distance method to use
    7. m: parameter for fuzzy c means clustering when calculating the weight matrix
    8. a_logit: log function parameter

    Output:
    1.w: weight matirx
    2. conitions: synthetic conditions based on n_conditions
    3. centroid_distance: distance from cells to centroids
    """
    np.random.seed(seed)  # set.seed
    conditions = [f"Condition{i}" for i in range(1, n_conditions + 1)]

    ## Find cluster center
    cluster_membership = adata.obs[pop_column]
    centroid_emb, cluster_names = find_centroid(X_emb, cluster_membership)

    centroid_distance = euclidean_distance(X_emb, centroid_emb)

    w = calculate_weights_centroid(centroid_distance, m)

    # Convert the NumPy array 'w' to a Pandas DataFrame
    w_df = pd.DataFrame(w)

    # Set column names from 'centroid_emb'
    w_df.columns = cluster_names

    # Set row names from 'X_emb'
    w_df.index = X_emb.index

    scaler = StandardScaler()
    w_scaled = scaler.fit_transform(w)
    w_logit = pd.DataFrame(
        log_function(w_scaled, a=a_logit), columns=w_df.columns, index=w_df.index
    )
    return w_logit, conditions
