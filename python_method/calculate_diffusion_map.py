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

def calculate_dm(adata,embedding_layer,n_dm):
    if adata.n_vars <= 50:
        adata.obsm[embedding_layer] = adata.X
    else:
        if embedding_layer not in adata.obsm:
            raise Exception("There is no embedding layer as input, please check your data.")
    
    dm_res = palantir.utils.run_diffusion_maps(adata, n_components=n_dm,pca_key = embedding_layer)
    adata_copy = adata.copy()
    adata_copy.obsm[embedding_layer] = adata.obsm["DM_EigenVectors"]
    
    return adata_copy