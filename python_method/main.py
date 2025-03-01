import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
from importlib import reload
import sys

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
import runMELD
import runMellon

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
    parser.add_argument('--mellon_d_method', type=str, help='The d_method for Mellon, can be "embedding" or "fractal"')
    parser.add_argument('--norm_density', type=str, help='when using d_method as fractal, do we want to use normalized density to calculate log fold change or use unnormalized')
    parser.add_argument('--hyperparameter', type=str, help='whether use computed hyperparameter for Mellon')
    parser.add_argument('--ls_factor', type=float, help='LS factor')
    parser.add_argument('--beta', type=float, help='Beta value')
    parser.add_argument('--k_meld', type=int, help='K MELD value')
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
    mellon_d_method = args.mellon_d_method
    norm_density = args.norm_density
    hyperparameter = args.hyperparameter
    ls_factor = args.ls_factor
    beta = args.beta
    k_meld = args.k_meld
    output_dir = args.output_dir
    #np.random.seed(seed)

    #pandas2ri.activate()

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
    
    output_dir = os.path.join(os.getcwd(), args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    adata = read_file.read_dataset(file_path, mode_embedding, pca_path, umap_path, n_components=10)


    if ds_type != "cluster":
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


            
        
        elif mode_select == "random":
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
            cell_type_df.to_csv(output_dir+'random_selected_cell_index'+mode_select+'.csv', index=False)

        
        performance_result_mellon = []
        
        performance_result_meld = []
        
        mellon_result_log = []
        mellon_result_zscore = []
        mellon_density_1 = []
        mellon_density_2 = []
        mellon_density_1_norm = []
        mellon_density_2_norm = []
        
        meld_sample_likelihood = []
        
        true_labels = []
        
        d_diff_cells = []
        mu_diff_cells = []
        ls_diff_cells = []
        cond1_prob = []
        cond2_prob = []
        
        print("start to work on iteration")
        
        if mode_select == "centroid":
            n_iteration = 1
        elif mode_select == "random":
            n_iteration = n_random_cell
        for i in range(n_iteration):
            print(i)
            print(hyperparameter)
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
            
            true_labels.append(list(adata.obs['true_labels']))
            cond1_prob.append(list(adata.obs["condition1_prob"]))
            cond2_prob.append(list(adata.obs["condition2_prob"]))

            
            if hyperparameter == "No":
            
                if n_dm <= 50:

                    condition1_dens,condition2_dens,condition1_dens_norm,condition2_dens_norm,log_fold_change_mean, zscores,d_list,ls_list,mu_list = runMellon.runMELLON(
                    adata,mellon_d_method, norm_density,"synth_labels" ,n_dm,ls_factor)

                elif n_dm >50:
                    condition1_dens,condition2_dens,condition1_dens_norm,condition2_dens_norm,log_fold_change_mean, zscores = runMellon.runMELLON(
                    adata, mellon_d_method, norm_density,"synth_labels" ,n_dm,n_dm,ls_factor)

            elif hyperparameter == "Yes":

                    condition1_dens,condition2_dens,condition1_dens_norm,condition2_dens_norm,log_fold_change_mean, zscores,d_list,ls_list,mu_list = runMellon.runMELLON_2(
                    adata, mellon_d_method, norm_density,"synth_labels" ,n_dm,ls_factor)

            mellon_thresholds = runMellon.threshold_mellon(zscores)
            print(mellon_thresholds)
            da_cell_mellon = runMellon.mellon2output(log_fold_change_mean, zscores, out_type="label", thresholds=mellon_thresholds)
            mellon_performance = runMellon.get_performance_df(adata,"true_labels",da_cell_mellon)
            mellon_performance.to_csv(output_dir+'mellon_performance_'+str(i)+'.csv', index=False)
            auc_mellon,prc_mellon = runMellon.evaluation(mellon_performance)

            performance_result_mellon.append(auc_mellon)
            performance_result_mellon.append(prc_mellon)

            mellon_result_log.append(log_fold_change_mean)
            mellon_result_zscore.append(zscores)
            mellon_density_1.append(condition1_dens)
            mellon_density_2.append(condition2_dens)
            mellon_density_1_norm.append(condition1_dens_norm)
            mellon_density_2_norm.append(condition2_dens_norm)
            
            d_diff_cells.append(d_list)
            mu_diff_cells.append(mu_list)
            ls_diff_cells.append(ls_list)
            
            
            sample_likelihoods_meld,samplem_meld = runMELD.runMELD(adata,k_meld,"synth_samples",'synth_labels', mode_embedding,beta)
            
            threshold_meld_i = runMELD.threshold_meld(adata,sample_likelihoods_meld)
            print(threshold_meld_i)
            
            da_cell_meld = runMELD.meld2output(sample_likelihoods_meld, out_type="label", thresholds=threshold_meld_i)
            meld_performance = runMELD.get_performance_df_meld(adata,"true_labels",da_cell_meld)
            meld_performance.to_csv(output_dir + 'meld_performance_'+str(i)+'.csv', index=False)
            
            
            auc_meld,prc_meld = runMELD.evaluation(meld_performance)
            
            performance_result_meld.append(auc_meld)
            performance_result_meld.append(prc_meld)
            
            meld_sample_likelihood.append(sample_likelihoods_meld)

    elif ds_type == "cluster":
        adata = cluster_dataset_synth_labels.add_synth_label_cluster_labels(adata,pop,seed, pop_enr,pop_column,n_conditions,n_batches, n_replicates,balance,cap_enr = None)
        adata = cluster_dataset_synth_labels.label_condition_and_rep_other(adata,n_replicates, n_batches,seed)
        if balance == "Yes":
            
            adata = cluster_dataset_synth_labels.quantile_assign_label(adata,pop_column,pop_enr,pop)
        elif balance == "No":
            adata = cluster_dataset_synth_labels.quantile_assign_label_old(adata,pop_column,pop_enr,pop)


        
        adata = cluster_dataset_synth_labels.add_batch_effect(adata, batch_col="synth_batches", norm_sd=batch_sd,seed = seed)
        true_labels = []
        cond1_prob = []
        cond2_prob = []
        true_labels.append(list(adata.obs['true_labels']))
        cond1_prob.append(list(adata.obs["condition1_prob"]))
        cond2_prob.append(list(adata.obs["condition2_prob"]))

        performance_result_mellon = []
        
        performance_result_meld = []
        
        mellon_result_log = []
        mellon_result_zscore = []
        mellon_density_1 = []
        mellon_density_2 = []
        mellon_density_1_norm = []
        mellon_density_2_norm = []
        
        meld_sample_likelihood = []
        
        
        
        d_diff_cells = []
        mu_diff_cells = []
        ls_diff_cells = []

        if mode_select == "centroid":
            n_iteration = 1
        elif mode_select == "random":
            n_iteration = n_random_cell
        for i in range(n_iteration):

            if hyperparameter == "No":
                
                    if n_dm <= 50:

                        condition1_dens,condition2_dens,condition1_dens_norm,condition2_dens_norm,log_fold_change_mean, zscores,d_list,ls_list,mu_list = runMellon.runMELLON(
                        adata,mellon_d_method, norm_density, "synth_labels" ,n_dm,ls_factor)

                    elif n_dm >50:
                        condition1_dens,condition2_dens,condition1_dens_norm,condition2_dens_norm,log_fold_change_mean, zscores = runMellon.runMELLON(
                        adata, mellon_d_method, norm_density,"synth_labels" ,n_dm,n_dm,ls_factor)

            elif hyperparameter == "Yes":

                        condition1_dens,condition2_dens,condition1_dens_norm,condition2_dens_norm,log_fold_change_mean, zscores,d_list,ls_list,mu_list = runMellon.runMELLON_2(
                        adata, mellon_d_method, norm_density,"synth_labels" ,n_dm,ls_factor)

            mellon_thresholds = runMellon.threshold_mellon(zscores)
            print(mellon_thresholds)
            da_cell_mellon = runMellon.mellon2output(log_fold_change_mean, zscores, out_type="label", thresholds=mellon_thresholds)
            mellon_performance = runMellon.get_performance_df(adata,"true_labels",da_cell_mellon)
            mellon_performance.to_csv(output_dir+'mellon_performance_'+str(i)+'.csv', index=False)
            auc_mellon,prc_mellon = runMellon.evaluation(mellon_performance)

            performance_result_mellon.append(auc_mellon)
            performance_result_mellon.append(prc_mellon)

            mellon_result_log.append(log_fold_change_mean)
            mellon_result_zscore.append(zscores)
            mellon_density_1.append(condition1_dens)
            mellon_density_2.append(condition2_dens)
            mellon_density_1_norm.append(condition1_dens_norm)
            mellon_density_2_norm.append(condition2_dens_norm)
                
            d_diff_cells.append(d_list)
            mu_diff_cells.append(mu_list)
            ls_diff_cells.append(ls_list)
                
                
            sample_likelihoods_meld,samplem_meld = runMELD.runMELD(adata,k_meld,"synth_samples",'synth_labels', mode_embedding,beta)
                
            threshold_meld_i = runMELD.threshold_meld(adata,sample_likelihoods_meld)
            print(threshold_meld_i)
                
            da_cell_meld = runMELD.meld2output(sample_likelihoods_meld, out_type="label", thresholds=threshold_meld_i)
            meld_performance = runMELD.get_performance_df_meld(adata,"true_labels",da_cell_meld)
            meld_performance.to_csv(output_dir + 'meld_performance_'+str(i)+'.csv', index=False)
                
                
            auc_meld,prc_meld = runMELD.evaluation(meld_performance)
                
            performance_result_meld.append(auc_meld)
            performance_result_meld.append(prc_meld)
                
            meld_sample_likelihood.append(sample_likelihoods_meld)

        
    
        # Ensure the output directory exists
       

    
    mellon_result_log_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(mellon_result_log)})
    mellon_result_zscore_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(mellon_result_zscore)})
    mellon_density_1_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(mellon_density_1)})
    mellon_density_2_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(mellon_density_2)})
    mellon_density_1_df_norm = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(mellon_density_1_norm)})
    mellon_density_2_df_norm = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(mellon_density_2_norm)})
    
    
    d_diff_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(d_diff_cells)})
    mu_diff_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(mu_diff_cells)})
    ls_diff_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(ls_diff_cells)})
    
    meld_result_sample_likelihood_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(meld_sample_likelihood)})
    true_labels_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(true_labels)})
    

    cond1_prob_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(cond1_prob)})
    cond2_prob_df = pd.DataFrame({f'cell_{i}': arr for i, arr in enumerate(cond2_prob)})
    
    mellon_result_log_df.to_csv(output_dir+"log_result_mellon_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    mellon_density_1_df.to_csv(output_dir+"mellon_density_condition1_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    mellon_density_2_df.to_csv(output_dir+"mellon_density_condition2_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    mellon_density_1_df_norm.to_csv(output_dir+"normalized_mellon_density_condition1_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    mellon_density_2_df_norm.to_csv(output_dir+"normalized_mellon_density_condition2_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    
    mellon_result_zscore_df.to_csv(output_dir+"zscore_result_mellon_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    meld_result_sample_likelihood_df.to_csv(output_dir+"meld_result_sample_likelihood_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    
    true_labels_df.to_csv(output_dir+"true_labels_"+str(n_dm)+"_"+str(n_iteration)+".txt")

    cond1_prob_df.to_csv(output_dir+"cond1_prob_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    cond2_prob_df.to_csv(output_dir+"cond2_prob_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    
    d_diff_df.to_csv(output_dir+"d_diff_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    
    ls_diff_df.to_csv(output_dir+"ls_diff_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    
    mu_diff_df.to_csv(output_dir+"mu_diff_"+str(n_dm)+"_"+str(n_iteration)+".txt")
   
        
    save_txt.save_list_to_file(performance_result_mellon, output_dir+"performance_result_mellon_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    save_txt.save_list_to_file(performance_result_meld, output_dir+"performance_result_meld_"+str(n_dm)+"_"+str(n_iteration)+".txt")
    
    #adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
    if 'X_pca_batch' in adata.obsm:
        del adata.obsm['X_pca_batch']
    adata.obs.columns = adata.obs.columns.astype(str)
    #adata.obsm.columns = adata.obsm.columns.astype(str)
    adata.var.columns = adata.var.columns.astype(str)
    #adata.write_h5ad(output_dir+"adata_"+str(n_dm)+"_"+str(n_iteration)+"_"+pop+".h5ad")

# Check if this script is being run directly
if __name__ == "__main__":
    main()