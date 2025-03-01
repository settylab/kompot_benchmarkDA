import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad


from sklearn.neighbors import NearestNeighbors

def find_random_cells_via_centroid(centroid_emb,X_emb,cluster_names,cluster_membership):
    
    """
    find the most nearest neighbor to the centroid to replace centroid with a real cell in the cluster
    Parameter:
    1. centroid_emb: centroid based on PCA or DiffusionMap, result fron find_centroid
    2. X_emb: cells based on PCA or DiffusionMap
    3. cluster_names: names for cell type represented by centroid_emb
    4. cluster_membership: adata.obs["pop_col"]
    
    Output: 
    1. centroid_nn_emb: nearest neighbor to the centroid
    """
    knn_temp = NearestNeighbors(n_neighbors=2)   ### This knn model is not for distance calculation, this knn model is just for finding nearest neighbor to centroid
    centroid_nn_emb = pd.DataFrame()
    
    centroid_emb_df = pd.DataFrame(centroid_emb, columns= X_emb.columns)
    centroid_index_list = cluster_names
    centroid_emb_df.index = centroid_index_list
    
    for cluster in cluster_names:
        cluster_points = X_emb[cluster_membership == cluster]
    
        subset_emb = pd.DataFrame(centroid_emb_df.loc[cluster,:]).T
        # assumes pca df as a dataFrame, can be changed later
        
        if isinstance(cluster_points, pd.DataFrame):
            cluster_points_df = cluster_points
        else:
            cluster_points_df = pd.DataFrame(cluster_points)
        
        
        combined_emb = pd.concat([subset_emb,cluster_points_df])
    
        knn_temp.fit(combined_emb)
        neighbor_distance,neighbors_cells = knn_temp.kneighbors(combined_emb)
        neighbor_distance,neighbors_cells = neighbor_distance[:, 1:],neighbors_cells[:, 1:]
        cell_index =neighbors_cells[0]
        cell_name = combined_emb.iloc[cell_index,:].index[0]
        cell_info = combined_emb.loc[cell_name,:].to_numpy()
        cell_info_df = pd.DataFrame(cell_info,columns = [cell_name])
        centroid_nn_emb = pd.concat([centroid_nn_emb, cell_info_df], axis=1)
    centroid_nn_emb = centroid_nn_emb.T
    centroid_nn_emb.columns = X_emb.columns
    
    return centroid_nn_emb



def random_sampling_larger(adata, centroid_nn_emb,pop_col,X_emb,cluster_membership,n_random_cells):
    
    """
    Parameters:
    1. adata: original data
    2. centroid_nn_emb: result from find_random_cells_via_centroid function
    3. pop_col: population column
    4. X_emb: cells on embedding space
    5. cluster_membership: adata.obs["pop_col"]
    
    Results:
    1. selected_cells: list of nearest neighbor cells of centroids
    2. random_selected_cells: dataframe of random cells and centroid nearest neighbors on embedding space
    3. flat_list_cell_name: list of all randomly selected_cells
    4. cell_type_dict: dictionary which include all random cells for each cell type
    """
    
    random_selected_cells = pd.DataFrame()
    cluster_names = []
    cell_index = []
    
    cell_type_dict = {}
    
    selected_cells = []
    for i in range(len(centroid_nn_emb.index)):
        cell_type_i = adata[adata.obs_names == centroid_nn_emb.index[i]].obs[pop_col][0]
        selected_cells.append(cell_type_i)
        
    for j in range(len(selected_cells)):
        cluster_points = X_emb[cluster_membership == selected_cells[j]]

        # Calculate the centroid and add it to the list
        random_cell_index = np.random.choice(len(cluster_points),size = n_random_cells-1,replace = False)
        
        random_cell = cluster_points.iloc[random_cell_index,:]
        
#         # add nearest neighbor of centroid to random cell list
        nn_cluster = pd.DataFrame(centroid_nn_emb.iloc[j,:]).T
        
        combined_cells = pd.concat([random_cell,nn_cluster])
        
        random_selected_cells = pd.concat([random_selected_cells,combined_cells])
        
        cell_name = list(cluster_points.iloc[random_cell_index,:].index)
        cell_name.append(str(centroid_nn_emb.iloc[j,:].name))

        cell_index.append(cell_name)
        cell_type_dict[selected_cells[j]] = cell_name
        
    flat_list_cell_name = []

#     # Loop through each list in 'nested_list' and add its items to 'flattened_list'
    for sublist in cell_index:
        for k in sublist:
            flat_list_cell_name.append(k)
            
    print(len(random_selected_cells))
        
    return selected_cells,random_selected_cells,flat_list_cell_name,cell_type_dict


def random_sampling_one(X_emb,cluster_membership):
    """
    When number of random cells equals to 1 and the nearest neighbor of centroid are not considered.
    Parameters: 
    1. X_emb: cells on embedding space
    2. cluster_membership: adata.obs[pop_col]
    
    Return:
    1. random_selected_cells: random selected cells, one cell for one cell type.
    2. cluster_names: names for each cell type represented by random_selected_cells
    3. cell_index: names for randomly selected cells
    """
    
    unique_clusters = np.unique(cluster_membership)
    random_selected_cells = []
    cluster_names = []
    cell_index = []
    
    cell_type_dict = {}

    for cluster in unique_clusters:

        cluster_points = X_emb[cluster_membership == cluster]

        # Calculate the centroid and add it to the list
        random_cell_index = np.random.choice(len(cluster_points))
        random_cell = np.array(cluster_points.iloc[random_cell_index,:])
        random_selected_cells.append(random_cell)
        cell_name = cluster_points.iloc[random_cell_index,:].name
        
        cluster_names.append(cluster)
        cell_index.append(cell_name)
        
        cell_type_dict[cluster] = [cell_name]
        
    random_selected_cells = np.array(random_selected_cells)
    
    
    return random_selected_cells,cluster_names,cell_index,cell_type_dict