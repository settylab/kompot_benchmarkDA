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
    Ensure that both PCA and DM embeddings are available with consistent batch correction.

    This function guarantees that:
    1. {layer_embedding}_batch exists (PCA with batch effects)
    2. If dm_comp > 0, DM_EigenVectors is computed FROM the batch-affected PCA
    3. All methods use the same batch-corrected embeddings

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

    batch_key = f"{layer_embedding}_batch"

    # 1. Ensure batch-corrected PCA exists
    if batch_key not in adata.obsm:
        logger.warning(f"{batch_key} not found in adata.obsm. This should have been created during preprocessing.")
        # Fallback: use clean embedding if batch-corrected doesn't exist
        if layer_embedding in adata.obsm:
            logger.warning(f"Using clean {layer_embedding} instead of batch-corrected version.")
            adata.obsm[batch_key] = adata.obsm[layer_embedding].copy()
        elif "X_pca" in adata.obsm:
            logger.warning(f"Using X_pca instead of {layer_embedding}.")
            adata.obsm[batch_key] = adata.obsm["X_pca"].copy()
        else:
            logger.error("No PCA embedding found. Computing PCA now.")
            sc.tl.pca(adata, n_comps=50)
            adata.obsm[batch_key] = adata.obsm["X_pca"].copy()

    # 2. Compute batch-corrected DM from batch-corrected PCA if needed
    dm_batch_key = "DM_EigenVectors_batch"
    if dm_comp > 0:
        # Check if batch-corrected DM exists with correct number of components
        if (dm_batch_key not in adata.obsm or
            adata.obsm[dm_batch_key].shape[1] != dm_comp):

            logger.info(f"Computing batch-corrected diffusion maps ({dm_comp} components) from batch-corrected PCA")

            # CRITICAL: Use batch-corrected PCA as input for DM computation
            # This will create DM_EigenVectors but we'll rename it to preserve the clean version
            temp_dm_key = "DM_EigenVectors"

            # Save the clean DM if it exists
            clean_dm = adata.obsm.get(temp_dm_key, None)

            palantir.utils.run_diffusion_maps(
                adata,
                n_components=dm_comp,
                pca_key=batch_key  # Use batch-corrected PCA!
            )

            # Move the batch-corrected DM to its proper location
            adata.obsm[dm_batch_key] = adata.obsm[temp_dm_key].copy()

            # Restore the clean DM if it existed
            if clean_dm is not None:
                adata.obsm[temp_dm_key] = clean_dm
            else:
                # If no clean DM existed, remove the temp one to avoid confusion
                del adata.obsm[temp_dm_key]

            logger.info(f"{dm_batch_key} shape: {adata.obsm[dm_batch_key].shape}")

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
        X = adata.obsm["DM_EigenVectors_batch"]
        embedding_name = "DM_EigenVectors_batch"
    else:
        X = adata.obsm["X_pca_batch"]
        embedding_name = "X_pca_batch"

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
        return "DM_EigenVectors_batch"
    else:
        return "X_pca_batch"