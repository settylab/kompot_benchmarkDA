import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad


#from scipy.spatial.distance import pdist, squareform, cdist

def calculate_weights_centroid(centroid_dist, m=2):
    
    """
    Parameters:
    1. centroid_dist: euclidean distance between cells and centroids
    2. m: parameter for fuzzy c means clustering
    
    return:
    1. w: weight matrix calculated based on centroid distance
    """
    
    n_rows, n_cols = centroid_dist.shape
    w = np.zeros_like(centroid_dist)
    
    for j in range(n_cols):
        for i in range(n_rows):
            w[i, j] = 1 / np.sum(centroid_dist[i, j] / centroid_dist[i, :]) ** (2 / (m - 1))
    
    return w


def calculate_weights_random_cell(random_cell_dist,m=2,eps = 1e-16):
    
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
    
    #if n_random_cell_per_pop > 1:
   #temp_dist = cdist(random_selected_cells, centroid_emb,metric = "euclidean")
    #print(temp_dist.shape)
    for j in range(n_cols):
        for i in range(n_rows):
            #if random_cell_dist[i,j] != 0:
                w[i, j] = 1 / np.sum((random_cell_dist[i, j] +eps) / (random_cell_dist[i, :] + eps)) ** (2 / (m - 1))
            #if random_cell_dist[i,j] == 0:
                #w[i, j] = 1 / np.sum(temp_dist[j, j//n_random_cell_per_pop] / centroid_dist[j, :]) ** (2 / (m - 1))
#     else:
#          for j in range(n_cols):
#             for i in range(n_rows):
#                 w[i, j] = 1 / np.sum(random_cell_dist[i, j] / centroid_dist[i, :]) ** (2 / (m - 1))
    
    return w