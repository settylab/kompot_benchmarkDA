import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import anndata as ad

import sys

import os

import re

import palantir

## Read anndata. For RDS, need to transfer it from SingleCellExperiment object 

def read_dataset(filepath,mode = "PCA", pca_file = None , umap_file = None, n_components = 10):
    """
    read adata
    Parameters:
    1. filepath : physical location of the anndata file
    2. pca_file: if the anndata file has separated PCA data which is not included in the anndata.obsm, read it. 
    3. mode: "PCA" or "DiffusionMap", indicate which embedding space to work on, default = "PCA"
    4. n_components: number of components for DiffusionMap
    
    The default for pca_file should be None.
    If the pca_file is None and "X_pca" not in adata.obsm, it will calculate PCA using scanpy.

    """
    adata = sc.read_h5ad(filepath)
    if mode == "PCA":
        if "X_pca" in adata.obsm:
            pca_data = adata.obsm["X_pca"]
        elif "X_pca" not in adata.obsm:
            if pca_file is not None:
                pca_data = pd.read_csv(pca_file,sep = "\t")
                adata.obsm["X_pca"] = pca_data.to_numpy()
                adata.uns['pca_var_names'] = pca_data.columns.tolist()
            else:
                sc.pp.normalize_total(adata)
                sc.pp.log1p(adata)
                #sc.pp.highly_variable_genes(adata,n_top_genes=2500)
                sc.pp.pca(adata,use_highly_variable=False, n_comps=50)
                sc.pp.pca(adata)
                
                
    elif mode == "DiffusionMap":
        if "X_pca" in adata.obsm:
            palantir.utils.run_diffusion_maps(adata,n_components)
        elif "X_pca" not in adata.obsm:
            if pca_file is not None:
                pca_data = pd.read_csv(pca_file,sep = "\t")
                adata.obsm["X_pca"] = pca_data.to_numpy()
                adata.uns['pca_var_names'] = pca_data.columns.tolist()
                palantir.utils.run_diffusion_maps(adata,n_components)
            else:
                sc.pp.normalize_total(adata)
                sc.pp.log1p(adata)
                #sc.pp.highly_variable_genes(adata,n_top_genes=2500)
                sc.pp.pca(adata,use_highly_variable=False, n_comps=50)
                sc.pp.pca(adata)
                palantir.utils.run_diffusion_maps(adata,n_components)
                
    if "X_umap" in adata.obsm:
        umap_data = adata.obsm["X_umap"]
    elif "X_umap" not in adata.obsm:
        if umap_file is not None:
            umap_data = pd.read_csv(umap_file,sep = "\t")
            adata.obsm["X_umap"] = umap_data.to_numpy()
            adata.uns['umap_var_names'] = umap_data.columns.tolist()
        else:
            sc.pp.neighbors(ad, use_rep="X_pca")
            sc.tl.umap(adata)
    
    return adata


## Get PCA dataframe for calculating distance.

def get_embedding_value(adata, embedding_name = "X_pca"):
    
    """
    Save X_pca in adata.obsm as an individual dataframe to use for synthetic labeling
    Parameters:
    1. adata: adata from read_dataset
    2. embedding_name: which obsm layer to work on, default is "X_pca", can be "DM_EigenVectors" or "X_pca"
    
    Return: emb_df: dataframe for pca data
    """
    emb_data = adata.obsm[embedding_name]
    emb_df = pd.DataFrame(emb_data, index=adata.obs_names)

    # Generate column names in the desired format (PC_1, PC_2, ...)
    emb_columns = [f"embedding_{i+1}" for i in range(emb_data.shape[1])]

    # Assign the generated column names to the DataFrame
    emb_df.columns = emb_columns
    
    return emb_df
