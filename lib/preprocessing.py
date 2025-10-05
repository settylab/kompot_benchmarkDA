import os
import warnings

import numpy as np
from scipy.stats import norm as normal


import anndata as ad
import scanpy as sc

import pandas as pd
from scipy.sparse import issparse

import scipy.sparse as sp


def replace_space_in_string(adata, pop_col):
    pop_names = adata.obs[pop_col].tolist()
    new_strings = []
    for s in pop_names:
        if " " in s:
            # print(f"Detected space in: '{s}'")
            s = s.replace(" ", "_")
        new_strings.append(s)
    adata.obs[pop_col] = new_strings
    adata.obs[pop_col] = adata.obs[pop_col].astype("category")
    return adata


def has_duplicate_rows(adata):

    if sp.issparse(adata.X):
        print("yes,sparse")
        dense_array = adata.X.toarray()
    else:
        print("not sparse, changing to dense array")
        dense_array = np.asarray(adata.X)
    array_of_arrays = dense_array
    seen = set()
    print("start to detecting duplications")
    for row in array_of_arrays:
        row_tuple = tuple(row)
        if row_tuple in seen:
            return True
        seen.add(row_tuple)
    return False


def remove_duplicate_cells(adata, keep="first"):
    """
    Remove duplicate cells from an AnnData object based on their expression matrix (adata.X).

    Duplicate cells are identified by comparing the rows of adata.X. Cells with identical
    expression profiles are considered duplicates, and only the first occurrence is kept
    (by default). The function returns a new AnnData object with duplicates removed.

    Parameters:
        adata (AnnData): The input AnnData object.
        keep (str): Which duplicate to keep. Options are:
                    - 'first' (default): Keep the first occurrence.
                    - 'last': Keep the last occurrence.
                    - False: Drop all duplicates.

    Returns:
        AnnData: A new AnnData object with duplicate cells removed.
    """

    # Convert the expression matrix to a DataFrame for easy duplicate detection.
    # Use adata.obs_names as the index and adata.var_names as the columns.
    if issparse(adata.X):
        X_df = pd.DataFrame(
            adata.X.toarray(), index=adata.obs_names, columns=adata.var_names
        )
    else:
        X_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)

    # Identify duplicated rows. The duplicated() method marks duplicates as True.
    duplicates = X_df.duplicated(keep=keep)
    num_duplicates = duplicates.sum()
    print(f"Number of duplicated cells: {num_duplicates}")

    # Keep only the unique rows (cells) based on expression
    unique_cell_names = X_df[~duplicates].index

    # Subset the AnnData object to only include unique cells.
    adata_unique = adata[unique_cell_names].copy()

    return adata_unique
