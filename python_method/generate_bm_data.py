import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad
from importlib import reload
import sys
from pathlib import Path
import os

from sklearn.preprocessing import StandardScaler


import condition_prob_centroid
import calculate_diffusion_map
import get_weight_matrix
import helper_functions
import read_file
# import runMELD
# import runMellon
import helper_functions

import argparse
import numpy as np


import sys
import os

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

    parser.add_argument('--n_conditions', type=int, help='Number of conditions')
    parser.add_argument('--n_replicates', type=int, help='Number of replicates')
    parser.add_argument('--n_batches', type=int, help='Number of batches')
    parser.add_argument('--seed', type=int, help='Seed for random number generation')
    parser.add_argument('--condition_balance', type=int, help='Condition balance')

    parser.add_argument('--m', type=float, help='M value')
    parser.add_argument('--a_logit', type=float, help='A_logit value')
    
    parser.add_argument('--mode_embedding', type=str, help='Embedding mode ,PCA or DiffusionMap')
    parser.add_argument('--layer_embedding', type=str, help='Layer embedding, X_pca or DM_EigenVectors')
    parser.add_argument('--n_dm', type=int, default=0, help='Number of diffusion components to compute on batch-simulated embeddings')
    parser.add_argument('--balance', type=str, help='whether we want to balance number of cells in each condition manually')
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
    n_conditions = args.n_conditions
    n_replicates = args.n_replicates
    n_batches = args.n_batches
    seed = args.seed
    condition_balance = args.condition_balance
    m = args.m
    a_logit = args.a_logit
    cap_enr = None  # remains unchanged
    mode_embedding = args.mode_embedding
    layer_embedding = args.layer_embedding
    n_dm = args.n_dm
    balance = args.balance
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
            f"Micromamba env path: {env_path}\n"
            "Please make sure you have R installed in the micromamba environment."
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
        # Use newer rpy2 API with context manager
        from rpy2.robjects.conversion import localconverter
        with localconverter(robjects.default_converter + pandas2ri.converter):
            return robjects.conversion.py2rpy(df)


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

    # Handle pandas2ri activation - suppress deprecation warning
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        try:
            pandas2ri.activate()
        except Exception:
            # If activation fails, continue without it
            pass
    

    adata = read_file.read_dataset(file_path,layer_embedding)
    output_dir = Path(output_dir)
    print(ds_type)
    print(output_dir)
    if ds_type != "cluster":
        print(ds_type)
        X_emb = read_file.get_embedding_value(adata,layer_embedding)
        print("X_emb")
        w_logit,conditions= get_weight_matrix.get_weight_matrix_centroid(adata, pop_column,seed, X_emb, n_conditions,m , a_logit)
            
        print("weight calculation done")
            
        enr_scores = condition_prob_centroid.generate_enr_prob(pop,pop_enr,w_logit)
            
        print("enr scores calculation done")
        enr_prob = condition_prob_centroid.normalize_enr_prob(w_logit,enr_scores, condition_balance)
            
        print("enr_prob calculation done")
        cond_probability = condition_prob_centroid.set_relevant_prob(enr_prob,pop_enr,pop, adata, pop_column)
        print("cond_probability done")
            
        w_logit.to_csv(output_dir / "weight_matrix.csv")
            
    ## Create a folder to save the future results
        
    

        print("Generating labels and batch-simulated embeddings")

        # Generate labels directly to output_dir (no iteration subdirectory)
        if not isinstance(cond_probability, pd.DataFrame):
            cond_probability_temp = pd.DataFrame(cond_probability,index = adata.obs_names)
        else:
            cond_probability_temp = cond_probability

        temp_cond = cond_probability_temp.iloc[:,0]

        conditions,cond_probability_df = synth_labels.cap_probabilities(adata,temp_cond,conditions, balance, cap_enr=None)


        adata = synth_labels.label_condition_and_rep_labels(adata,cond_probability_df,seed)


        adata = synth_labels.label_condition_and_rep_other(adata,n_replicates, n_batches,seed)

        if balance == "Yes":

            adata = synth_labels.quantile_assign_label(adata,pop_column,pop_enr,pop)
        elif balance == "No":
            adata = synth_labels.quantile_assign_label_old(adata,pop_column,pop_enr,pop)

        if mode_embedding == "PCA":
            adata = synth_labels.add_batch_effect_pca(adata,layer_embedding, batch_col="synth_batches", norm_sd=batch_sd,seed = seed)

            # Compute DM on batch-simulated embeddings if needed
            if n_dm > 0:
                print(f"Computing DM ({n_dm} components) on batch-simulated embeddings...")
                adata = calculate_diffusion_map.calculate_dm(adata, f"{layer_embedding}_batch", n_dm)
                print("DM computation on batch-simulated embeddings done")

            X_pca = pd.DataFrame(adata.obsm[f"{layer_embedding}_batch"],index=adata.obs_names)
            print("done")
        else:
            adata = synth_labels.add_batch_effect_pca(adata,layer_embedding, batch_col="synth_batches", norm_sd=batch_sd,seed = seed)

            # Compute DM on batch-simulated embeddings if needed
            if n_dm > 0:
                print(f"Computing DM ({n_dm} components) on batch-simulated embeddings...")
                adata = calculate_diffusion_map.calculate_dm(adata, f"{layer_embedding}_batch", n_dm)
                print("DM computation on batch-simulated embeddings done")

            X_pca = pd.DataFrame(adata.obsm[f"{layer_embedding}_batch"],index=adata.obs_names)


        # Assuming `adata` is your AnnData object
        obs_df = adata.obs.copy()

        # Add the index as a column named "rowname"
        obs_df['rowname'] = obs_df.index

        cols = ['rowname'] + [col for col in obs_df.columns if col != 'rowname']
        obs_df = obs_df[cols]

        obs_df.to_csv(output_dir / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}.coldata.csv',index = False)

        str_batch = str(batch_sd)
        int_batch = helper_functions.convert_number_str(str_batch)

        # Save X_pca_batch (for PCA-based methods like MELD, CNA, R methods)
        X_pca["rowname"] = obs_df.index
        cols = ['rowname'] + [col for col in X_pca.columns if col != 'rowname']
        X_pca = X_pca[cols]
        X_pca.to_csv(output_dir / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.emb.csv', index = False,float_format='%.16f')

        # Save DM_EigenVectors_batch separately (for DM-based methods like Mellon, Kompot)
        if n_dm > 0:
            X_dm = pd.DataFrame(adata.obsm["DM_EigenVectors_batch"], index=adata.obs_names)
            X_dm["rowname"] = obs_df.index
            cols = ['rowname'] + [col for col in X_dm.columns if col != 'rowname']
            X_dm = X_dm[cols]
            X_dm.to_csv(output_dir / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.emb.dm.csv', index = False,float_format='%.16f')

    
    #"cluster" or
    elif ds_type ==  "cluster":
        adata = cluster_dataset_synth_labels.add_synth_label_cluster_labels(adata,pop,seed, pop_enr,pop_column,n_conditions,balance,cap_enr = None)
        adata = cluster_dataset_synth_labels.label_condition_and_rep_other(adata,n_replicates, n_batches,seed)
        if balance == "Yes":

            adata = cluster_dataset_synth_labels.quantile_assign_label(adata,pop_column,pop_enr,pop)
        elif balance == "No":
            adata = cluster_dataset_synth_labels.quantile_assign_label_old(adata,pop_column,pop_enr,pop)

        if mode_embedding == "PCA":
            adata = cluster_dataset_synth_labels.add_batch_effect_pca(adata, layer_embedding,batch_col="synth_batches", norm_sd=batch_sd,seed = seed)

            # Compute DM on batch-simulated embeddings if needed
            if n_dm > 0:
                print(f"Computing DM ({n_dm} components) on batch-simulated embeddings...")
                adata = calculate_diffusion_map.calculate_dm(adata, f"{layer_embedding}_batch", n_dm)
                print("DM computation on batch-simulated embeddings done")

            X_pca = pd.DataFrame(adata.obsm[f"{layer_embedding}_batch"],index=adata.obs_names)
        else:
            adata = cluster_dataset_synth_labels.add_batch_effect_pca(adata,layer_embedding, batch_col="synth_batches", norm_sd=batch_sd,seed = seed)

            # Compute DM on batch-simulated embeddings if needed
            if n_dm > 0:
                print(f"Computing DM ({n_dm} components) on batch-simulated embeddings...")
                adata = calculate_diffusion_map.calculate_dm(adata, f"{layer_embedding}_batch", n_dm)
                print("DM computation on batch-simulated embeddings done")

            X_pca = pd.DataFrame(adata.obsm[f"{layer_embedding}_batch"],index=adata.obs_names)

        obs_df = adata.obs.copy()

        # Add the index as a column named "rowname"
        obs_df['rowname'] = obs_df.index
        cols = ['rowname'] + [col for col in obs_df.columns if col != 'rowname']
        obs_df = obs_df[cols]
        obs_df.to_csv(output_dir / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}.coldata.csv',index =False)
        str_batch = str(batch_sd)
        int_batch = helper_functions.convert_number_str(str_batch)

        # Save X_pca_batch (for PCA-based methods like MELD, CNA, R methods)
        X_pca["rowname"] = obs_df.index
        cols = ['rowname'] + [col for col in X_pca.columns if col != 'rowname']
        X_pca = X_pca[cols]
        X_pca.to_csv(output_dir / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.emb.csv', index = False,float_format='%.16f')

        # Save DM_EigenVectors_batch separately (for DM-based methods like Mellon, Kompot)
        if n_dm > 0:
            X_dm = pd.DataFrame(adata.obsm["DM_EigenVectors_batch"], index=adata.obs_names)
            X_dm["rowname"] = obs_df.index
            cols = ['rowname'] + [col for col in X_dm.columns if col != 'rowname']
            X_dm = X_dm[cols]
            X_dm.to_csv(output_dir / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}.emb.dm.csv', index = False,float_format='%.16f')

        
        


if __name__ == "__main__":
    main()