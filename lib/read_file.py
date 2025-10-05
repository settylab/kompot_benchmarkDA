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


def read_dataset(filepath, layer_embedding):
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
    if adata.n_vars <= 50:
        if "X_x" in adata.obsm:
            adata.obsm[layer_embedding] = adata.obsm["X_x"]
        else:
            if not isinstance(adata.X, np.ndarray):
                adata.obsm[layer_embedding] = adata.X.toarray()
            else:
                adata.obsm[layer_embedding] = adata.X

    return adata


## Get PCA dataframe for calculating distance.


def get_embedding_value(adata, embedding_name="X_pca"):
    """
    Save X_pca in adata.obsm as an individual dataframe to use for synthetic labeling
    Parameters:
    1. adata: adata from read_dataset
    2. embedding_name: which obsm layer to work on, default is "X_pca", can be "DM_EigenVectors" or "X_pca"

    Return: emb_df: dataframe for pca data
    """

    if embedding_name not in adata.obsm:
        raise ValueError(
            f"Embedding '{embedding_name}' not found in adata.obsm. "
            "Please make sure it has been calculated."
        )

    emb_data = adata.obsm[embedding_name]
    emb_df = pd.DataFrame(emb_data, index=adata.obs_names)

    # Generate column names in the desired format (embedding_1, embedding_2, ...)
    emb_df.columns = [f"embedding_{i+1}" for i in range(emb_data.shape[1])]

    return emb_df
