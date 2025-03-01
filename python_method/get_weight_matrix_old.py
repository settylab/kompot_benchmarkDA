import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.preprocessing import StandardScaler

from find_centroid import find_centroid
from euclidean_distance import euclidean_distance, cosine_distance
from weight_calculation_old import calculate_weights_centroid,calculate_weights_random_cell
from helper_functions import normalization,log_function

from random_selected import find_random_cells_via_centroid, random_sampling_larger, random_sampling_one

from knn_distance import knn_distance_calculation


def get_weight_matrix_random_old(adata, pop,pop_column,seed, X_emb, mode = "euclidean", n_random_cells = 10,n_conditions = 2,knn_k = 15,m = 2, a_logit = 0.5):
    
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
    3. random_cell_distance: distance from cells to centroids
    4. cell_type_dict: dictionary for 
    """
    
    np.random.seed(seed) # set.seed
    conditions = [f"Condition{i}" for i in range(1, n_conditions + 1)]
    
    ## Find cluster center
    cluster_membership = adata.obs[pop_column]
    centroid_emb,cluster_names = find_centroid(X_emb, cluster_membership)
    
    centroid_distance = euclidean_distance(X_emb, centroid_emb)
    
    if n_random_cells > 1:
        centroid_nn_emb = find_random_cells_via_centroid(centroid_emb,X_emb,cluster_names,cluster_membership)
        selected_cells,samples_emb, samples_name,cell_type_dict = random_sampling_larger(adata, centroid_nn_emb,pop_column,X_emb,cluster_membership,n_random_cells)
        if mode == "euclidean":
            distance = euclidean_distance(X_emb, samples_emb)
        elif mode == "KNN":
            distance = knn_distance_calculation(samples_emb,cluster_names,X_emb,knn_k,samples_name,"random")
            
        random_cell_distance = pd.DataFrame(distance)
        w = calculate_weights_random_cell(distance,samples_emb,centroid_emb, centroid_distance,n_random_cells,m=2)
        
            
        # Convert the NumPy array 'w' to a Pandas DataFrame
        w_df = pd.DataFrame(w)

        # Set column names from 'centroid_emb'
        w_df.columns = samples_name
        # Set row names from 'X_emb'
        w_df.index = X_emb.index

    elif n_random_cells  ==  1:
        random_selected_cells,cluster_names,cell_index,cell_type_dict = random_sampling_one(X_emb,cluster_membership)
        if mode == "euclidean":
            distance = euclidean_distance(X_emb, random_selected_cells)
        elif mode == "KNN":
            distance = knn_distance_calculation(random_selected_cells,cluster_names,X_emb,knn_k,cell_index,"random")

        random_cell_distance = pd.DataFrame(distance)
        w = calculate_weights_random_cell(distance,random_selected_cells,centroid_emb, centroid_distance,n_random_cells,m=2)
    
        # Convert the NumPy array 'w' to a Pandas DataFrame
        w_df = pd.DataFrame(w)

        # Set column names from 'centroid_emb'
        w_df.columns = cell_index
        # Set row names from 'X_emb'
        w_df.index = X_emb.index

    w_scaled = normalization(w)
    
    w_logit = pd.DataFrame(log_function(w_scaled, a=a_logit), columns=w_df.columns, index=w_df.index)

    random_cell_distance = pd.DataFrame(random_cell_distance,index = w_df.index, columns = w_df.columns)
    
    return w_logit,conditions,random_cell_distance,cell_type_dict
    
    
    
    
    
    ### If we choose to work on centroid
def get_weight_matrix_centroid_old(adata, pop,pop_column,seed, X_emb, n_conditions, knn_k,mode = "euclidean",m = 2, a_logit = 0.5,cell_index = None):
    
    
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
    np.random.seed(seed) # set.seed
    conditions = [f"Condition{i}" for i in range(1, n_conditions + 1)]
    
    ## Find cluster center
    cluster_membership = adata.obs[pop_column]
    centroid_emb,cluster_names = find_centroid(X_emb, cluster_membership)
    
    if mode == "euclidean":
        centroid_distance = euclidean_distance(X_emb, centroid_emb)
    elif mode == "KNN":
        centroid_distance = knn_distance_calculation(centroid_emb,cluster_names,X_emb,knn_k,None,mode == "centroid")
    elif mode == "cosine":
        centroid_distance = cosine_distance(centroid_emb,cluster_names,X_emb,knn_k,None,mode == "centroid")

    w = calculate_weights_centroid(centroid_distance,m)
    
    # Convert the NumPy array 'w' to a Pandas DataFrame
    w_df = pd.DataFrame(w)

    # Set column names from 'centroid_emb'
    w_df.columns = cluster_names

    # Set row names from 'X_emb'
    w_df.index = X_emb.index
    
    scaler = StandardScaler()
    w_scaled = scaler.fit_transform(w)
    w_logit = pd.DataFrame(log_function(w_scaled, a=a_logit), columns=w_df.columns, index=w_df.index)
    return w_logit,conditions,centroid_distance