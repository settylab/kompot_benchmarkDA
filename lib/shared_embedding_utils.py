"""
Shared embedding utilities to ensure consistency across all DA methods.
This module provides standardized functions for computing embeddings with batch effects.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import palantir
import logging

logger = logging.getLogger(__name__)

def ensure_batch_corrected_embeddings(adata, layer_embedding="X_pca", dm_comp=10):
    """
    Ensure that both PCA and DM embeddings are available.

    Since data_loader loads batch-affected embeddings directly into standard keys
    (X_pca, DM_EigenVectors), this function now just validates they exist and
    optionally computes DM if not present.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    layer_embedding : str, default "X_pca"
        Base embedding layer name
    dm_comp : int, default 10
        Number of diffusion map components (0 = PCA only)

    Returns:
    --------
    adata : AnnData
        Modified AnnData with consistent embeddings
    """

    # 1. Ensure PCA exists
    if layer_embedding not in adata.obsm:
        logger.warning(f"{layer_embedding} not found in adata.obsm.")
        # Fallback: use X_pca if a different layer was requested
        if layer_embedding != "X_pca" and "X_pca" in adata.obsm:
            logger.warning(f"Using X_pca instead of {layer_embedding}.")
            # Don't copy - just let methods fall back to X_pca
        else:
            logger.error("No PCA embedding found. Computing PCA now.")
            sc.tl.pca(adata, n_comps=50)

    # 2. Compute DM from PCA if needed
    dm_key = "DM_EigenVectors"
    if dm_comp > 0:
        # Check if DM exists with correct number of components
        if (dm_key not in adata.obsm or
            adata.obsm[dm_key].shape[1] != dm_comp):

            logger.info(f"Computing diffusion maps ({dm_comp} components) from PCA")

            # Compute DM from the PCA (which already has batch effects if loaded from CSV)
            palantir.utils.run_diffusion_maps(
                adata,
                n_components=dm_comp,
                pca_key=layer_embedding  # Use the PCA embedding (has batch effects)
            )

            logger.info(f"{dm_key} shape: {adata.obsm[dm_key].shape}")

    return adata

def get_embedding_for_method(adata, use_dm=True, dm_comp=10):
    """
    Get the appropriate embedding matrix for a DA method.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    use_dm : bool, default True
        Whether to use diffusion maps (True) or PCA (False)
    dm_comp : int, default 10
        Number of DM components (only used if use_dm=True)

    Returns:
    --------
    X : np.ndarray
        Embedding matrix to use for the method
    embedding_name : str
        Name of the embedding used
    """

    # Ensure consistent embeddings exist
    adata = ensure_batch_corrected_embeddings(adata, dm_comp=dm_comp)

    if use_dm and dm_comp > 0:
        X = adata.obsm["DM_EigenVectors"]
        embedding_name = "DM_EigenVectors"
    else:
        X = adata.obsm["X_pca"]
        embedding_name = "X_pca"

    # Convert to numpy if needed
    if not isinstance(X, np.ndarray):
        X = np.array(X)

    logger.info(f"Using embedding: {embedding_name}, shape: {X.shape}")

    return X, embedding_name

def get_obsm_key_for_method(use_dm=True):
    """
    Get the obsm key that should be used for a DA method.

    Parameters:
    -----------
    use_dm : bool, default True
        Whether to use diffusion maps (True) or PCA (False)

    Returns:
    --------
    str
        The obsm key to use
    """
    if use_dm:
        return "DM_EigenVectors"
    else:
        return "X_pca"