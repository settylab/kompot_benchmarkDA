import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad

import networkx as nx

from sklearn.neighbors import NearestNeighbors

def knn_distance_calculation(representative_cell_emb,cluster_names, X_emb, knn_k = 15,cell_index = None,mode = "centroid"):
    
    """
    Calculation knn_distance using dijsktra algorithm
    
    Parameter:
    1. representative_cell_emb: randomly selected cells or centroids on pca or diffusionmap
    2. cluster_names: cell names for each centroid
    3. mode: "centroid" or "random"
    4. cell_index: cell names for randomly selected real cells
    5. X_emb: cells on pca or diffusionmap
    6. knn_k: a paramter for K Nearest Neighbor, number of nearest neighbors.
    
    Output:
    1. knn_distance_new: numpy array of numpy arrays for knn distance

    
    """
    
    k=knn_k
    ## Create knn model
    knn = NearestNeighbors(n_neighbors=k+1) # To avoid the function recognize the first nearest neighbor as itself
    combined_emb = pd.DataFrame()
    # assumes pca df as a dataFrame, can be changed later
    representative_cell_emb_df = pd.DataFrame(representative_cell_emb, columns= X_emb.columns)
    if mode == "centroid":
        representative_cell_emb_df.index = cluster_names
        combined_emb = pd.concat([representative_cell_emb_df,X_emb])
        final_cell_index = combined_emb.index
    elif mode == "random":
        representative_cell_emb_df.index = cell_index
        emb_temp = X_emb.drop(cell_index, axis=0)
        combined_emb = pd.concat([representative_cell_emb_df,emb_temp])
    
        final_cell_index = combined_emb.index
    
    # Combined_pca at here is like the reordered dataframe
    
    knn.fit(combined_emb)
    neighbor_distance,neighbors_cells = knn.kneighbors(combined_emb)
    
    neighbor_distance,neighbors_cells = neighbor_distance[:, 1:],neighbors_cells[:, 1:]
    
    G = nx.Graph()
    for i in range(len(list(combined_emb.index))):
        for j in neighbors_cells[i]:
            if i!= j:
                G.add_edge(i,j,weight = neighbor_distance[i][np.where(neighbors_cells[i] == j)[0][0]])
            
    shortest_distance = []
    for j in range(len(cell_index)):   
        shorted_distance_j = nx.single_source_dijkstra_path_length(G, source=j)
        shortest_distance.append(shorted_distance_j)

    if mode == "centroid":
        shortest_distance_filtered = []  # Remove centrioid to centroid distance
        for i in shortest_distance:
            new_i = {k: v for k, v in i.items() if k >= len(centroid_emb)}
            sorted_dict = {k: new_i[k] for k in sorted(new_i)}
            shortest_distance_filtered.append(sorted_dict)
    
        df_list = [pd.DataFrame.from_dict(d, orient='index', columns=[f'Column{i}']) for i, d in enumerate(shortest_distance_filtered, start=1)]
        knn_distance = pd.concat(df_list, axis=1)
    elif mode == "random":
        shortest_distance_filtered = []  # Remove centrioid to centroid distance
        for i in shortest_distance:
            sorted_dict = {k: i[k] for k in sorted(i)}
            shortest_distance_filtered.append(sorted_dict)
    
        df_list = [pd.DataFrame.from_dict(d, orient='index', columns=[f'Column{i}']) for i, d in enumerate(shortest_distance_filtered, start=1)]
        knn_distance = pd.concat(df_list, axis=1)
  
    # # if the graph is not 
    # if len(df_list) != len(combined_emb):
    #     raise ValueError("The graph G is not fully connected, please consider increasing k value")
    
    
    knn_distance.index = final_cell_index
    #knn_distance.columns = cluster_names
    new_knn_distance_df = knn_distance.loc[list(X_emb.index)]
     
    knn_distance_new = new_knn_distance_df.to_numpy()
    return  knn_distance_new