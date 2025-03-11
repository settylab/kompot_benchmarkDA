import read_file

import runMELD



import argparse
import numpy as np


import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
import sys
import helper_functions
from pathlib import Path

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

    parser.add_argument('--beta', type=float, help='Beta value')
    parser.add_argument('--k_meld', type=int, help='K MELD value')

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
    layer_embedding = args.layer_embedding
    input_file = args.input_file
    package = args.package

    #output_filepath = args.output_filepath
    beta = args.beta
    k_meld = args.k_meld

    output_dir = args.output_dir


    output_dir = Path(output_dir)
    input_file = Path(input_file)
    str_batch = str(batch_sd)
    int_batch = helper_functions.convert_number_str(str_batch)

    
    for i in range(1):
        iteration_directory = input_file / f'iteration_{i}'
        output_dir_i = output_dir / f'iteration_{i}'
        adata = read_file.read_dataset(file_path,layer_embedding)

        adata.obsm[f"{layer_embedding}_batch"] = np.array(pd.read_csv(iteration_directory/f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.emb.csv', index_col=0))

        adata.obs = pd.read_csv(iteration_directory/f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}.coldata.csv', index_col=0)
        
        sample_likelihoods_meld,samplem = runMELD.runMELD(adata,k_meld,"synth_samples",'synth_labels', layer_embedding,beta)

        df_meld = pd.DataFrame(sample_likelihoods_meld, columns=[f"col_{i}" for i in range(sample_likelihoods_meld.reshape(-1,1).shape[1])],index=adata.obs_names)
        
        df_meld.to_csv(output_dir_i / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}_package_performance.DAresults.{package}.csv')
        # meld_auroc_score,meld_auroc_score_neg, meld_auroc_score_pos = runMELD.modified_auroc(adata,sample_likelihoods_meld)
        # meld_auprc_score,meld_auprc_score_neg, meld_auprc_score_pos = runMELD.modified_auprc(adata,sample_likelihoods_meld)
        # meld_performance_dict = {"auroc":[meld_auroc_score],"auroc_neg":[meld_auroc_score_neg],"auroc_pos":[meld_auroc_score_pos],"auprc":[meld_auprc_score],"auprc_neg":[meld_auprc_score_neg],"auprc_pos":[meld_auprc_score_pos]}
        # meld_performance = pd.DataFrame(meld_performance_dict)
        # threshold_meld_i = runMELD.threshold_meld(adata,sample_likelihoods_meld)
        # print(threshold_meld_i)
            
        # da_cell_meld = runMELD.meld2output(sample_likelihoods_meld, out_type="label", thresholds=threshold_meld_i)
        # meld_performance = runMELD.get_performance_df_meld(adata,"true_labels",da_cell_meld,threshold_meld_i)

        #meld_performance.to_csv(output_dir_i / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.DAresults.{package}.modified.csv', index=False)

if __name__ == "__main__":
    main()