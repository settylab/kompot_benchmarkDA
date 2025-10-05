from scipy.spatial.distance import pdist, squareform, cdist
import numpy as np
import pandas as pd


def euclidean_distance(X_emb, centroid_emb):
    """
    X_emb: result of get_embedding_value function, can be on "PCA" space or "DiffusionMap" space
    centroid_emb: centroids on embedding space

    Return:
    c_distance: The distance from cells in X_emb to centroids. The row number should be number of cells and the column number should be number of centroids
    """
    c_distance = cdist(X_emb, centroid_emb, metric="euclidean")
    return c_distance
