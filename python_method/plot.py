import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad

def clip_distance(distance_matrix,selected_columns,weight_matrix,adata):
    adata_temp = adata.copy()
    distance_matrix_df = pd.DataFrame(distance_matrix,index = weight_matrix.index,columns = weight_matrix.columns)
    distance_matrix_temp = distance_matrix_df[selected_columns]
    weight_matrix_temp = weight_matrix[selected_columns]
    
    distance_matrix_log = np.log10(distance_matrix_temp)
    distance_matrix_clipped = np.clip(
    distance_matrix_log, *np.quantile(distance_matrix_log, [0.01, 1]))
    distance_df = pd.DataFrame(distance_matrix_clipped,index = weight_matrix_temp.index,columns = weight_matrix_temp.columns)
    adata_temp.obs = pd.concat([adata_temp.obs,distance_df],axis = 1)
    return adata_temp


def distance_plot(adata,selected_columns):
    sc.pl.embedding(adata, "umap", color=selected_columns,color_map = "Spectral_r")
    adata.obs.drop(selected_columns,axis = 1,inplace = True)
    
def weight_plot(adata, weight, selected_columns):
    w_cell = weight[selected_columns]
    adata.obs = pd.concat([adata.obs,w_cell],axis = 1)
    sc.pl.embedding(adata, "umap", color=selected_columns,color_map = "Spectral_r")
    adata.obs.drop(selected_columns,axis = 1,inplace = True)
    
import seaborn as sns
def box_plot(list_of_performance, y_axis):
    # Creating the box plots
    plt.boxplot(list_of_performance)

    # Generating x-axis labels dynamically based on the number of lists
    labels = ['List {}'.format(i+1) for i in range(len(list_of_performance))]

    # Setting the x-axis labels
    plt.xticks(range(1, len(list_of_performance)+1), labels)

    # Adding title and labels for clarity
    plt.title('Comparison of performance scores across different diffusion maps')
    plt.ylabel(y_axis)
    

