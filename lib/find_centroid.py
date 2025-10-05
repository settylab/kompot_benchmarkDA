import numpy as np
import anndata


def find_centroid(X_emb, cluster_membership):
    """
    find centroid for each type in pop_col
    Parameters:
    1. X_emb: PCA or DiffusionMap
    2. cluster_membership: which cluster that each cell belongs to, need to be extracted from anndata.obs[pop_col]

    Return:
    1. centroid_emb: centroids based on embedding space
    2. cluster_names: names for each embedding
    """

    unique_clusters = np.unique(cluster_membership)
    centroid_emb = []
    cluster_names = []

    for cluster in unique_clusters:
        cluster_points = X_emb[cluster_membership == cluster]

        # Calculate the centroid and add it to the list
        cluster_centroid = cluster_points.mean(axis=0)
        centroid_emb.append(cluster_centroid)

        cluster_names.append(cluster)
    centroid_emb = np.array(centroid_emb)

    return centroid_emb, cluster_names
