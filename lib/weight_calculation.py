import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from lib.constants import FUZZY_CMEANS_M, FUZZY_CMEANS_EPS

# from scipy.spatial.distance import pdist, squareform, cdist


def calculate_weights_centroid(centroid_dist, m=FUZZY_CMEANS_M, eps=FUZZY_CMEANS_EPS):
    """
    Calculate fuzzy membership weights using fuzzy c-means algorithm.

    The weight w[i,j] represents the fuzzy membership of cell i to centroid j.
    Formula: w[i,j] = 1 / sum_k[(d[i,j] / d[i,k])^(2/(m-1))]

    Parameters:
    -----------
    centroid_dist : array-like, shape (n_cells, n_centroids)
        Euclidean distance matrix between cells and centroids
    m : float, default=2
        Fuzziness parameter for fuzzy c-means clustering.
        m > 1 required. Larger m → fuzzier memberships.
        m=2 is standard choice (moderate fuzziness).
    eps : float, default=1e-10
        Small value added to distances to prevent division by zero
        when distances are very small or cells are exactly at centroids

    Returns:
    --------
    w : array, shape (n_cells, n_centroids)
        Fuzzy membership weight matrix. Each row sums to ~1.0.
        w[i,j] closer to 1 means cell i belongs more to centroid j.

    Notes:
    ------
    - Edge case: If all distances for a cell are equal, returns uniform weights
    - Numerical stability: eps added to prevent division by zero
    - Input validation: m must be > 1 for valid fuzzy c-means
    """
    if m <= 1:
        raise ValueError(f"Fuzziness parameter m must be > 1, got m={m}")

    n_rows, n_cols = centroid_dist.shape
    w = np.zeros_like(centroid_dist, dtype=float)

    # Add eps for numerical stability
    centroid_dist_stable = centroid_dist + eps

    exponent = 2 / (m - 1)

    for j in range(n_cols):
        for i in range(n_rows):
            # Calculate ratio of distances
            distance_ratios = centroid_dist_stable[i, j] / centroid_dist_stable[i, :]

            # Sum of ratios raised to exponent
            ratio_sum = np.sum(distance_ratios ** exponent)

            # Handle edge case: if all distances equal, weights should be uniform
            if ratio_sum == 0 or not np.isfinite(ratio_sum):
                w[i, j] = 1.0 / n_cols
            else:
                w[i, j] = 1.0 / ratio_sum

    return w


def calculate_weights_random_cell(random_cell_dist, m=2, eps=1e-16):
    """
    If we select to use real cells instead of centroid, the distance between cell and cell itself will be 0, which will cause error in the fuzze c means clustering.
    If replacing zero with a very smallest value closed to 0, the weight will be very large, which will influenece the normalization and sigmoid function result.
    At here, I tried to replace 0 with its distance to the centroid.

    If the n_random_cell_per_pop is larger than 1, it means that we have 1 nearest neighbor cell and at least 1 real randomly selected cells.
    If the n_random_cell_per_pop equals to 1, it means that we only have 1 randomly selected cells.

    Parameters:
    1. random_cell_dist: the distance from cells to random cells, calculated by euclidean distance or knn distance
    2. random_selected_cells: the randomly selected cells on embedding
    3. centroid_emb: centroids on embedding space
    4. entroid_dist: cells to centroids distances
    5. n_random_cell_per_pop: number of randommly selected cells + 1 nearest neighbor to the centroid
    6. eps: small value close to zero to avoid error

    i is for each cell
    j is for each cell type

    Output:
    1. w: weight matrix

    Others:
    temp_dist: distance from random slected cells and their distance to the centroid of that cell type
    """

    n_rows, n_cols = random_cell_dist.shape
    w = np.zeros_like(random_cell_dist)

    # if n_random_cell_per_pop > 1:
    # temp_dist = cdist(random_selected_cells, centroid_emb,metric = "euclidean")
    # print(temp_dist.shape)
    for j in range(n_cols):
        for i in range(n_rows):
            # if random_cell_dist[i,j] != 0:
            w[i, j] = 1 / np.sum(
                (random_cell_dist[i, j] + eps) / (random_cell_dist[i, :] + eps)
            ) ** (2 / (m - 1))
        # if random_cell_dist[i,j] == 0:
        # w[i, j] = 1 / np.sum(temp_dist[j, j//n_random_cell_per_pop] / centroid_dist[j, :]) ** (2 / (m - 1))
    #     else:
    #          for j in range(n_cols):
    #             for i in range(n_rows):
    #                 w[i, j] = 1 / np.sum(random_cell_dist[i, j] / centroid_dist[i, :]) ** (2 / (m - 1))

    return w
