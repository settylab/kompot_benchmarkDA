import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
from sklearn.preprocessing import StandardScaler

from find_centroid import find_centroid
from euclidean_distance import euclidean_distance
from weight_calculation import calculate_weights_centroid
from helper_functions import log_function



## If we choose to work on random cells
# def get_weight_matrix_random(adata, pop,pop_column,seed, X_emb, mode = "euclidean", n_random_cells = 10,n_conditions = 2,knn_k = 15,m = 2, a_logit = 0.5):
    
#     """
#     Calculated final weight matrix when the distance is calculated from cells to centroids
#     1. adata: original adata to work on
#     2. pop: which population
#     3. pop_column: the columns for population in adata.obs
#     4. seed: random seed
#     5. X_emb: cells on embedding space
#     6. mode: can be "euclidean" or "KNN", which distance method to use
#     7. m: parameter for fuzzy c means clustering when calculating the weight matrix
#     8. a_logit: log function parameter
    
#     Output:
#     1.w: weight matirx
#     2. conitions: synthetic conditions based on n_conditions
#     3. random_cell_distance: distance from cells to centroids
#     4. cell_type_dict: dictionary for 
#     """
    
#     np.random.seed(seed) # set.seed
#     conditions = [f"Condition{i}" for i in range(1, n_conditions + 1)]
    
#     ## Find cluster center
#     cluster_membership = adata.obs[pop_column]
#     centroid_emb,cluster_names = find_centroid(X_emb, cluster_membership)
    
#     centroid_distance = euclidean_distance(X_emb, centroid_emb)
    
#     if n_random_cells > 1:
#         centroid_nn_emb = find_random_cells_via_centroid(centroid_emb,X_emb,cluster_names,cluster_membership)
#         selected_cells,samples_emb, samples_name,cell_type_dict = random_sampling_larger(adata, centroid_nn_emb,pop_column,X_emb,cluster_membership,n_random_cells)
#         print(len(samples_name))
#         if mode == "euclidean":
#             distance = euclidean_distance(X_emb, samples_emb)
#         elif mode == "KNN":
#             distance = knn_distance_calculation(samples_emb,cluster_names,X_emb,knn_k,samples_name,"random")
#         if mode == "cosine":
#             distance = cosine_distance(X_emb, samples_emb)
#         random_cell_distance = pd.DataFrame(distance)
        
#         w_logit = pd.DataFrame(index = X_emb.index)
#         print(len(list(set(list(adata.obs[pop_column])))))
#         print(distance.shape)
#         print(np.transpose(distance).shape)
#         distance_t = np.transpose(distance)
#         for i in range(n_random_cells):
#             # Placeholder for the selected column data
#             selected_distance = []
#             sample_names_j = []
#             # Loop through each group of 10 arrays
#             for j in range(len(list(set(list(adata.obs[pop_column]))))):
#                 # For each group, select the column from the first array
#                 selected_column = distance_t[j*n_random_cells+i]
#                 sample_name_ij = samples_name[j*n_random_cells+i]
#                 # Append the selected column data to the results list
#                 selected_distance.append(selected_column)
#                 sample_names_j.append(sample_name_ij)
#             # Convert selected_columns to a NumPy array if needed
#             selected_distance_array = np.transpose(np.array(selected_distance))
#             print(selected_distance_array.shape)
#             w_temp = calculate_weights_random_cell(selected_distance_array,m=2,eps = 1e-16)
#             print(w_temp.shape)
#                     # Convert the NumPy array 'w' to a Pandas DataFrame
#             w_df_temp = pd.DataFrame(w_temp)
#             w_df_temp.columns = sample_names_j
#             w_df_temp.index = X_emb.index
            
#             w_scaled_temp = normalization(w_temp)
    
#             w_logit_temp = pd.DataFrame(log_function(w_scaled_temp, a=a_logit), columns=w_df_temp.columns, index=w_df_temp.index)
        
            
#             w_logit = pd.concat([w_logit,w_logit_temp],axis = 1)
#         w_df = pd.DataFrame(columns = samples_name,index = X_emb.index)


#     elif n_random_cells  ==  1:
#         random_selected_cells,cluster_names,cell_index,cell_type_dict = random_sampling_one(X_emb,cluster_membership)
#         if mode == "euclidean":
#             distance = euclidean_distance(X_emb, random_selected_cells)
#         elif mode == "KNN":
#             distance = knn_distance_calculation(random_selected_cells,cluster_names,X_emb,knn_k,cell_index,"random")
#         if mode == "cosine":
#             distance = cosine_distance(X_emb, random_selected_cells)
            
#         random_cell_distance = pd.DataFrame(distance)
#         w = calculate_weights_random_cell(distance,m=2,eps = 1e-16)
    
#         # Convert the NumPy array 'w' to a Pandas DataFrame
#         w_df = pd.DataFrame(w)

#         # Set column names from 'centroid_emb'
#         w_df.columns = cell_index
#         # Set row names from 'X_emb'
#         w_df.index = X_emb.index

#         w_scaled = normalization(w)

#         w_logit = pd.DataFrame(log_function(w_scaled, a=a_logit), columns=w_df.columns, index=w_df.index)

#     random_cell_distance = pd.DataFrame(random_cell_distance,index = w_df.index, columns = w_df.columns)
    
#     return w_logit,conditions,random_cell_distance,cell_type_dict
    
    
    
    
    
    
    ### If we choose to work on centroid
def get_weight_matrix_centroid(adata, pop_column,seed, X_emb, n_conditions ,m = 2, a_logit = 0.5):
    
    
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
    
    centroid_distance = euclidean_distance(X_emb, centroid_emb)

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
    return w_logit,conditions