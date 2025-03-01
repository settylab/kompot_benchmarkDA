import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
from importlib import reload
import sys
from pathlib import Path
import os

import re

from scipy.spatial.distance import pdist, squareform, cdist
from sklearn.preprocessing import StandardScaler


import condition_prob_centroid,condition_prob_random
import euclidean_distance
import find_centroid
import get_weight_matrix
import helper_functions
import knn_distance
import plot
import random_selected
import read_file
# import runMELD
# import runMellon
import helper_functions

import weight_calculation
import save_txt

import argparse
import numpy as np


import sys
import os

def main():
        # Create the parser
    parser = argparse.ArgumentParser(description='Your script description.')

    # Add arguments
    parser.add_argument('--file_path', type=str, help='Path to the file')
    parser.add_argument('--mode_distance', type=str, help='Distance mode, euclidean or KNN')
    parser.add_argument('--pca_path', type=str, help='PCA path')
    parser.add_argument('--umap_path', type=str, help='UMAP path')
    parser.add_argument('--pop', type=str, help='Population')
    parser.add_argument('--pop_enr', type=float, help='Population enrichment')
    parser.add_argument('--pop_column', type=str, help='Population column')
    parser.add_argument('--mode_select', type=str, help='Mode select, centroid or random')

    parser.add_argument('--ds_type', type=str, help='type of dataset')
    parser.add_argument('--batch_sd', type=float, help='batch standard deviation')

    parser.add_argument('--n_conditions', type=int, help='Number of conditions')
    parser.add_argument('--n_replicates', type=int, help='Number of replicates')
    parser.add_argument('--n_batches', type=int, help='Number of batches')
    parser.add_argument('--seed', type=int, help='Seed for random number generation')
    parser.add_argument('--condition_balance', type=int, help='Condition balance')

    parser.add_argument('--m', type=float, help='M value')
    parser.add_argument('--a_logit', type=float, help='A_logit value')
    parser.add_argument('--n_random_cell', type=int, help='Number of random cells')
    
    parser.add_argument('--mode_embedding', type=str, help='Embedding mode ,PCA or DiffusionMap')
    parser.add_argument('--knn_k', type=int, help='KNN K value')
    parser.add_argument('--layer_embedding', type=str, help='Layer embedding, X_pca or DM_EigenVectors')
    parser.add_argument('--n_dm', type=int, help='Number of diffusion component for Mellon value')
    parser.add_argument('--balance', type=str, help='whether we want to balance number of cells in each condition manually')
    #parser.add_argument('--mellon_d_method', type=str, help='The d_method for Mellon, can be "embedding" or "fractal"')
    # parser.add_argument('--norm_density', type=str, help='when using d_method as fractal, do we want to use normalized density to calculate log fold change or use unnormalized')
    # parser.add_argument('--hyperparameter', type=str, help='whether use computed hyperparameter for Mellon')
    # parser.add_argument('--ls_factor', type=float, help='LS factor')
    # parser.add_argument('--beta', type=float, help='Beta value')
    # parser.add_argument('--k_meld', type=int, help='K MELD value')
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')

    # Parse arguments
    args = parser.parse_args()

    # Accessing arguments (example)
    file_path = args.file_path
    mode_distance = args.mode_distance
    pca_path = args.pca_path
    umap_path = args.umap_path
    pop = args.pop
    pop_enr = args.pop_enr
    pop_column = args.pop_column
    mode_select = args.mode_select
    ds_type = args.ds_type
    batch_sd = args.batch_sd
    n_conditions = args.n_conditions
    n_replicates = args.n_replicates
    n_batches = args.n_batches
    seed = args.seed
    condition_balance = args.condition_balance
    m = args.m
    a_logit = args.a_logit
    n_random_cell = args.n_random_cell
    cap_enr = None  # remains unchanged
    #output_filepath = args.output_filepath
    mode_embedding = args.mode_embedding
    knn_k = args.knn_k
    layer_embedding = args.layer_embedding
    n_dm = args.n_dm
    balance = args.balance
    #mellon_d_method = args.mellon_d_method
    # norm_density = args.norm_density
    # hyperparameter = args.hyperparameter
    # ls_factor = args.ls_factor
    # beta = args.beta
    # k_meld = args.k_meld
    output_dir = args.output_dir



    if "/app/software/R" in os.environ["PATH"]:
        print(
        "R is loaded in the environment. Please make sure this is not part of "
        "your jupyter server start script and restart it from a clean command "
        "line environment."
    )
    else:
            print("R does not seem to be loaded and you are good to go.")
        


        # Make sure no R module is loaded
    if "/app/software/R" in os.environ["PATH"]:
           raise Exception("An R module seems to be loaded.")

        # Get the path to the python executable
    python_executable_path = sys.executable

        # Extract the path to the environment from the path to the python executable
    env_path = os.path.dirname(os.path.dirname(python_executable_path))

    print(
            f"Conda env path: {env_path}\n"
            "Please make sure you have R installed in the conda environment."
        )
    print(env_path)
    os.environ['R_HOME'] = os.path.join(env_path, 'lib', 'R')
    
    import rpy2.robjects as robjects
    from rpy2.robjects import pandas2ri, r, ListVector

    #import rpy2.robjects as ro
    from rpy2.robjects.packages import importr
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

        # Convert Pandas DataFrame to R DataFrame
    def pandas_to_r_dataframe(df):
        pandas2ri.activate()
        return pandas2ri.py2rpy(df)


   # Function to get the library paths used by R
    def get_r_lib_paths():
        r_command = """
        .libPaths()
        """
        return robjects.r(r_command)

    # Execute the function and store the result
    r_lib_paths = get_r_lib_paths()

    # Print the first library path used by R (assuming there's at least one path)
    if r_lib_paths is not None:
        print("R is using packages from " + r_lib_paths[0])
    else:
        print("No library paths found.")
    
    import synth_labels
    import cluster_dataset_synth_labels

    pandas2ri.activate()
    

    adata = read_file.read_dataset(file_path, mode_embedding, pca_path, umap_path, n_components=10)
    output_dir = Path(output_dir)
    if ds_type != "cluster_balanced":
        X_emb = read_file.get_embedding_value(adata,layer_embedding)
        
        if mode_select == "centroid":
            if mode_distance == "euclidean":
                w_logit,conditions,centroid_distance = get_weight_matrix.get_weight_matrix_centroid(adata, pop,pop_column,seed, X_emb, n_conditions, knn_k, mode_distance,m , a_logit,cell_index = None)
            elif mode_distance == "KNN":
                w_logit,conditions,centroid_distance = get_weight_matrix.get_weight_matrix_centroid(adata, pop,pop_column,seed, X_emb, n_conditions, knn_k, mode_distance,m , a_logit,cell_index = None)
            elif mode_distance == "cosine":
                w_logit,conditions,centroid_distance = get_weight_matrix.get_weight_matrix_centroid(adata, pop,pop_column,seed, X_emb, n_conditions, knn_k, mode_distance,m , a_logit,cell_index = None)
            
            print("weight calculation done")
            
            enr_scores = condition_prob_centroid.generate_enr_prob(pop,pop_enr,w_logit)
            
            print("enr scores calculation done")
            enr_prob = condition_prob_centroid.normalize_enr_prob(w_logit,enr_scores, condition_balance)
            
            print("enr_prob calculation done")
            cond_probability = condition_prob_centroid.set_relevant_prob(enr_prob,pop_enr,pop, adata, pop_column)
            print("cond_probability done")
            
            w_logit.to_csv(output_dir / "weight_matrix.csv")
            
        
        elif mode_select == "random":
            if n_random_cell >= len(list(adata[adata.obs[pop_column] == pop].obs_names)):
                n_random_cell == len(list(adata[adata.obs[pop_column] == pop].obs_names))
            if mode_distance == "euclidean":
                w_logit,conditions,random_cell_distance,cell_type_dict = get_weight_matrix.get_weight_matrix_random(adata, pop,pop_column,seed, X_emb, mode_distance, n_random_cell,n_conditions,knn_k,m, a_logit)
            elif mode_distance == "KNN":
                w_logit,conditions,random_cell_distance,cell_type_dict = get_weight_matrix.get_weight_matrix_random(adata, pop,pop_column,seed, X_emb, mode_distance, n_random_cell,n_conditions,knn_k,m, a_logit)
            elif mode_distance == "cosine":
                w_logit,conditions,random_cell_distance,cell_type_dict = get_weight_matrix.get_weight_matrix_random(adata, pop,pop_column,seed, X_emb, mode_distance, n_random_cell,n_conditions,knn_k,m, a_logit)
                
            print("weight calculation done")
            enr_scores = condition_prob_random.generate_enr_prob(pop,pop_enr,w_logit,cell_type_dict)
            print("enr scores calculation done")
            enr_prob = condition_prob_random.normalize_enr_prob(w_logit,enr_scores, condition_balance)
            print("enr_prob calculation done")
            cond_probability = condition_prob_random.set_relevant_prob(enr_prob,pop_enr,pop, adata, pop_column,cell_type_dict)
            print("cond_probability done")
            cell_type_df = pd.DataFrame(dict([(k, pd.Series(v)) for k,v in cell_type_dict.items()]))
            cell_type_df.to_csv(output_dir / f'random_selected_cell_index_{mode_select}.csv', index=False)

            w_logit.to_csv(output_dir / "weight_matrix.csv")
    ## Create a folder to save the future results
        
    

        print("start to work on iteration")

        n_iteration = 0
        print(mode_select)
        if mode_select == "centroid":
                n_iteration = n_iteration + 1
        elif mode_select == "random":
                n_iteration = n_iteration + n_random_cell
        print(n_iteration)
        for i in range(n_iteration):
            print(i)
            iteration_directory = output_dir / f'iteration_{i}'
            iteration_directory.mkdir(parents=True, exist_ok=True)
            if not isinstance(cond_probability, pd.DataFrame):
        
                cond_probability_temp = pd.DataFrame(cond_probability,index = adata.obs_names)
            else:
                cond_probability_temp = cond_probability
                
            temp_cond = cond_probability_temp.iloc[:,i]
                
            conditions,cond_probability_df = synth_labels.cap_probabilities(adata,temp_cond,conditions, balance, cap_enr=None)
                
                
            adata = synth_labels.label_condition_and_rep_labels(adata,cond_probability_df,conditions,n_replicates, n_batches,seed)


            adata = synth_labels.label_condition_and_rep_other(adata,conditions,n_replicates, n_batches,seed)
                
                #adata = synth_labels.add_new_synthetic_data(adata, synth_label,synth_samples,synth_batches_df,cond_probability_df)
                
            if balance == "Yes":
                
                adata = synth_labels.quantile_assign_label(adata,pop_column,pop_enr,pop)
            elif balance == "No":
                adata = synth_labels.quantile_assign_label_old(adata,pop_column,pop_enr,pop)

            adata = synth_labels.add_batch_effect(adata, batch_col="synth_batches", norm_sd=batch_sd,seed = seed)
            X_pca = pd.DataFrame(adata.obsm["X_pca_batch"],index=adata.obs_names)
            

            # Assuming `adata` is your AnnData object
            obs_df = adata.obs.copy()

            # Add the index as a column named "rowname"
            obs_df['rowname'] = obs_df.index

            cols = ['rowname'] + [col for col in obs_df.columns if col != 'rowname']
            obs_df = obs_df[cols]

            obs_df.to_csv(iteration_directory / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}.coldata.csv',index = False)

            str_batch = str(batch_sd)
            int_batch = helper_functions.convert_number_str(str_batch)

            X_pca["rowname"] = obs_df.index
            cols = ['rowname'] + [col for col in X_pca.columns if col != 'rowname']
            X_pca = X_pca[cols]
            X_pca.to_csv(iteration_directory / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.pca.csv', index = False,float_format='%.16f')

            if mode_embedding == "DiffusionMap":
                X_dm = pd.DataFrame(adata.obsm["DM_EigenVectors"],index=adata.obs_names)
                X_dm["rowname"] = obs_df.index
                cols = ['rowname'] + [col for col in X_dm.columns if col != 'rowname']
                X_dm = X_dm[cols]
                X_dm.to_csv(iteration_directory / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.DM.csv', index = False,float_format='%.16f')
    
    
    #"cluster" or
    elif ds_type ==  "cluster_balanced":
        i = 0
        iteration_directory = output_dir / f'iteration_{i}'
        iteration_directory.mkdir(parents=True, exist_ok=True)
        adata = cluster_dataset_synth_labels.add_synth_label_cluster_labels(adata,pop,seed, pop_enr,pop_column,n_conditions,n_batches, n_replicates,balance,cap_enr = None)
        adata = cluster_dataset_synth_labels.label_condition_and_rep_other(adata,n_replicates, n_batches,seed)
        if balance == "Yes":
            
            adata = cluster_dataset_synth_labels.quantile_assign_label(adata,pop_column,pop_enr,pop)
        elif balance == "No":
            adata = cluster_dataset_synth_labels.quantile_assign_label_old(adata,pop_column,pop_enr,pop)
        adata = cluster_dataset_synth_labels.add_batch_effect(adata, batch_col="synth_batches", norm_sd=batch_sd,seed = seed)
        X_pca = pd.DataFrame(adata.obsm["X_pca_batch"],index=adata.obs_names)


        obs_df = adata.obs.copy()

            # Add the index as a column named "rowname"
        obs_df['rowname'] = obs_df.index
        cols = ['rowname'] + [col for col in obs_df.columns if col != 'rowname']
        obs_df = obs_df[cols]
        obs_df.to_csv(iteration_directory / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}.coldata.csv',index =False)
        str_batch = str(batch_sd)
        int_batch = helper_functions.convert_number_str(str_batch)
        X_pca["rowname"] = obs_df.index
        cols = ['rowname'] + [col for col in X_pca.columns if col != 'rowname']
        X_pca = X_pca[cols]
        X_pca.to_csv(iteration_directory / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.pca.csv', index = False,float_format='%.16f')

        if mode_embedding == "DiffusionMap":
            X_dm = pd.DataFrame(adata.obsm["DM_EigenVectors"],index=adata.obs_names)
            X_dm["rowname"] = obs_df.index
            cols = ['rowname'] + [col for col in X_dm.columns if col != 'rowname']
            X_dm = X_dm[cols]
            X_dm.to_csv(iteration_directory / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.DM.csv', index = False,float_format='%.16f')

        
        


if __name__ == "__main__":
    main()