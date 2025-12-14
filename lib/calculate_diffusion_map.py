import os
import warnings

import numpy as np
#from scipy.stats import norm as normal

#import matplotlib
#import matplotlib.pyplot as plt
#from matplotlib.lines import Line2D

import anndata as ad
import scanpy as sc
import palantir

def calculate_dm(adata, embedding_layer, n_dm):
    """
    Compute diffusion map on the specified embedding layer.

    Parameters:
    - adata: AnnData object
    - embedding_layer: Name of embedding to use as input (e.g., 'X_pca', 'X_pca_batch')
    - n_dm: Number of diffusion components to compute

    Returns:
    - adata with DM_EigenVectors (or DM_EigenVectors_batch) added to obsm
    - Original embedding layer is PRESERVED
    """
    # For CyTOF data (n_vars <= 50), ensure embedding_layer exists
    if adata.n_vars <= 50:
        if embedding_layer not in adata.obsm:
            if not isinstance(adata.X, np.ndarray):
                adata.obsm[embedding_layer] = adata.X.toarray()
            else:
                adata.obsm[embedding_layer] = adata.X
    else:
        if embedding_layer not in adata.obsm:
            raise Exception(f"Embedding layer '{embedding_layer}' not found in adata.obsm. Available: {list(adata.obsm.keys())}")

    # Determine output keys based on input embedding
    # If computing DM on batch-simulated embeddings, add _batch suffix to all outputs
    if embedding_layer.endswith('_batch'):
        kernel_key = 'DM_Kernel_batch'
        sim_key = 'DM_Similarity_batch'
        eigval_key = 'DM_EigenValues_batch'
        eigvec_key = 'DM_EigenVectors_batch'
    else:
        kernel_key = 'DM_Kernel'
        sim_key = 'DM_Similarity'
        eigval_key = 'DM_EigenValues'
        eigvec_key = 'DM_EigenVectors'

    # Run diffusion maps with custom output keys
    dm_res = palantir.utils.run_diffusion_maps(
        adata,
        n_components=n_dm,
        pca_key=embedding_layer,
        kernel_key=kernel_key,
        sim_key=sim_key,
        eigval_key=eigval_key,
        eigvec_key=eigvec_key
    )

    return adata