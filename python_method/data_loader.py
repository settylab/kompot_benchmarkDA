"""
Unified data loading module for all differential abundance methods.
Provides consistent data loading patterns across all methods.
"""

import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path
import helper_functions

def load_dataset(file_path, label_directory, ds_type, pop, pop_enr, seed, batch_sd, layer_embedding, use_dm=False):
    """
    Load a dataset and associated metadata for differential abundance analysis.

    Parameters:
    -----------
    file_path : str
        Path to the main dataset file
    label_directory : str or Path
        Directory containing label files
    ds_type : str
        Dataset type identifier
    pop : str
        Population identifier
    pop_enr : float
        Population enrichment value
    seed : int
        Random seed value
    batch_sd : float
        Batch standard deviation
    layer_embedding : str
        Layer embedding name (e.g., 'X_pca', 'DM_EigenVectors')
    use_dm : bool, optional
        If True, load DM embeddings (.emb.dm.csv), otherwise load PCA embeddings (.emb.csv)
        Default is False (load PCA embeddings)

    Returns:
    --------
    adata : AnnData
        Annotated data object with loaded metadata and embeddings
    """
    # Convert to Path objects
    label_directory = Path(label_directory)

    # Load the main dataset
    adata = read_dataset(file_path)

    # Convert batch_sd to string format
    str_batch = str(batch_sd)
    int_batch = helper_functions.convert_number_str(str_batch)

    # Determine which embedding file to load
    if use_dm:
        # Load DM embeddings for DM-based methods (Mellon, Kompot)
        emb_suffix = '.emb.dm.csv'
        embedding_key = 'DM_EigenVectors_batch'
    else:
        # Load PCA embeddings for PCA-based methods (MELD, CNA, R methods)
        emb_suffix = '.emb.csv'
        embedding_key = f'{layer_embedding}_batch'

    # Load embedding data
    adata.obsm[embedding_key] = np.array(
        pd.read_csv(
            label_directory / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}{emb_suffix}',
            index_col=0
        )
    )

    # Load observation metadata
    adata.obs = pd.read_csv(
        label_directory / f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}.coldata.csv',
        index_col=0
    )

    return adata

def read_dataset(filename):
    """
    Read a dataset file using the appropriate format.
    
    Parameters:
    -----------
    filename : str
        Path to the dataset file
    
    Returns:
    --------
    adata : AnnData
        Annotated data object
    """
    # Use scanpy's read function which handles multiple formats
    return sc.read_h5ad(filename)

def save_results(result_df, output_dir, ds_type, pop, pop_enr, seed, batch_sd, package, suffix=""):
    """
    Save results to a CSV file.
    
    Parameters:
    -----------
    result_df : DataFrame
        Results to save
    output_dir : str or Path
        Directory to save results to
    ds_type : str
        Dataset type identifier
    pop : str
        Population identifier
    pop_enr : float
        Population enrichment value
    seed : int
        Random seed value
    batch_sd : float
        Batch standard deviation
    package : str
        Package/method identifier
    suffix : str, optional
        Additional suffix for the filename
    """
    # Convert to Path object
    output_dir = Path(output_dir)
    
    # Convert batch_sd to string format
    str_batch = str(batch_sd)
    int_batch = helper_functions.convert_number_str(str_batch)
    
    # Create the output filename
    filename = f'benchmark_{ds_type}_pop_{pop}_enr{pop_enr}_seed{seed}_batchEffect{int_batch}{suffix}.DAresults.{package}.csv'
    
    # Save to CSV
    result_df.to_csv(output_dir / filename)