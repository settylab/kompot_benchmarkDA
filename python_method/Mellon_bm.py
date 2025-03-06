import read_file

import runMellon

from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score

import argparse
import numpy as np


import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import sys
from pathlib import Path
import helper_functions

import os

import re

def main():
        # Create the parser
    parser = argparse.ArgumentParser(description='Your script description.')

    # Add arguments
    parser.add_argument('--file_path', type=str, help='Path to the file')
    parser.add_argument('--pop', type=str, help='Population')
    parser.add_argument('--pop_enr', type=float, help='Population enrichment')
    parser.add_argument('--pop_column', type=str, help='Population column')

    parser.add_argument('--ds_type', type=str, help='type of dataset')
    parser.add_argument('--batch_sd', type=float, help='batch standard deviation')
    parser.add_argument('--input_file', type=str, help='input_pca_file_path')
    parser.add_argument('--package', type=str, help='which package is using')


    parser.add_argument('--seed', type=int, help='Seed for random number generation')

    parser.add_argument('--layer_embedding', type=str, help='Layer embedding, X_pca or DM_EigenVectors')
    parser.add_argument('--n_dm', type=int, help='Number of diffusion component for Mellon value')

    parser.add_argument('--mellon_d_method', type=str, help='The d_method for Mellon, can be "embedding" or "fractal"')
    parser.add_argument('--norm_density', type=str, help='when using d_method as fractal, do we want to use normalized density to calculate log fold change or use unnormalized')
    parser.add_argument('--hyperparameter', type=str, help='whether use computed hyperparameter for Mellon')
    parser.add_argument('--corrected', type=str, help='whether correct the log fold change mean of Mellon')
    parser.add_argument('--ls_factor', type=float, help='LS factor')
    parser.add_argument('--ls_mode', type=str, help='LS mode, if the embedding is on PCA, ls_mode = "PCA", if the embedding is on DM, the mode is on DM')


    parser.add_argument('--output_dir', type=str, required=True, help='Output directory path')

    # Parse arguments
    args = parser.parse_args()

    # Accessing arguments (example)
    file_path = args.file_path

    pop = args.pop
    pop_enr = args.pop_enr
    pop_column = args.pop_column

    ds_type = args.ds_type
    batch_sd = args.batch_sd

    seed = args.seed

    input_file = args.input_file
    package = args.package

    #output_filepath = args.output_filepath
    layer_embedding = args.layer_embedding
    n_dm = args.n_dm

    mellon_d_method = args.mellon_d_method
    norm_density = args.norm_density
    hyperparameter = args.hyperparameter
    corrected = args.corrected
    ls_factor = args.ls_factor
    ls_mode = args.ls_mode

    output_dir = args.output_dir


    output_dir = Path(output_dir)
    input_file = Path(input_file)
    str_batch = str(batch_sd)
    int_batch = helper_functions.convert_number_str(str_batch)
    
    #adata_orig = read_file.read_dataset(file_path, mode_embedding, pca_path, umap_path, n_components=10)


    for i in range(1):
        iteration_directory = input_file / f'iteration_{i}'
        output_dir_i = output_dir / f'iteration_{i}'
        adata = read_file.read_dataset(file_path)
        adata.obsm[f"{layer_embedding}_batch"] = np.array(pd.read_csv(iteration_directory/f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.emb.csv', index_col=0))

        adata.obs = pd.read_csv(iteration_directory/f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}.coldata.csv', index_col=0)
        if hyperparameter == "No":
    
            log_fold_change_mean, zscores = runMellon.runMELLON(
                    adata, mellon_d_method, norm_density,"synth_labels" ,n_dm,ls_factor)

        elif hyperparameter == "Yes":

            log_fold_change_mean, zscores = runMellon.runMELLON_synchronized(
                        adata, mellon_d_method, norm_density,corrected,"synth_labels" ,n_dm,ls_factor,ls_mode)
        
        

        df_mellon_lfc = pd.DataFrame(log_fold_change_mean, columns=[f"col_{i}" for i in range(log_fold_change_mean.reshape(-1,1).shape[1])],index=adata.obs_names)
        df_mellon_zscore = pd.DataFrame(zscores, columns=[f"col_{i}" for i in range(zscores.reshape(-1,1).shape[1])],index=adata.obs_names)
        # mellon_auroc_score,mellon_auroc_score_neg, mellon_auroc_score_pos = runMellon.modified_auroc(adata,log_fold_change_mean)
        # mellon_auprc_score,mellon_auprc_score_neg, mellon_auprc_score_pos = runMellon.modified_auprc(adata,log_fold_change_mean)
        # mellon_performance_dict = {"auroc":[mellon_auroc_score],"auroc_neg":[mellon_auroc_score_neg],"auroc_pos":[mellon_auroc_score_pos],"auprc":[mellon_auprc_score],"auprc_neg":[mellon_auprc_score_neg],"auprc_pos":[mellon_auprc_score_pos]}
        # mellon_performance = pd.DataFrame(mellon_performance_dict)
        # mellon_thresholds = runMellon.threshold_mellon(zscores)
        # print(mellon_thresholds)
        # da_cell_mellon = runMellon.mellon2output(log_fold_change_mean, zscores, out_type="label", thresholds=mellon_thresholds)
        # mellon_performance = runMellon.get_performance_df(adata,"true_labels",da_cell_mellon,mellon_thresholds)
        #mellon_performance.to_csv(output_dir_i / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.DAresults.{package}.modified.csv', index=False)

        df_mellon_lfc.to_csv(output_dir_i / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}_package_performance.DAresults.{package}.csv')
        df_mellon_zscore.to_csv(output_dir_i / "mellon_zscore.csv")
if __name__ == "__main__":
    main()