import os
import warnings

import numpy as np
from scipy.stats import norm as normal

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import anndata
import scanpy as sc
import palantir
import mellon

import pandas as pd

import sys
import warnings
import anndata
import numpy as np
import scanpy as sc
import os.path as osp

import cna
from multianndata import MultiAnnData
import shared_embedding_utils

def runCNA_func(
    adata: anndata.AnnData,
    k: int,
    sample_col: str,
    label_col: str,
    encode_label_dict: dict = {'Condition1': 0, "Condition2": 1},
    batch_col: str = None,
    layer_embedding: str = "X_pca",
    use_dm: bool = False,
    dm_comp: int = 10,
):
    # Ensure consistent embeddings
    adata = shared_embedding_utils.ensure_batch_corrected_embeddings(
        adata, layer_embedding=layer_embedding, dm_comp=dm_comp
    )

    # Get the appropriate embedding matrix for building kNN graph
    X, embedding_name = shared_embedding_utils.get_embedding_for_method(
        adata, use_dm=use_dm, dm_comp=dm_comp
    )

    # Store the embedding matrix for scanpy's neighbor computation
    adata.obsm['X_cna_embedding'] = X

    # Build the kNN graph using the standardized embedding
    sc.pp.neighbors(adata, n_neighbors=k, use_rep='X_cna_embedding')

    # create multi-anndata and convert condition and batch labels to numeric vars
    adata.obs['sample_id'] = adata.obs[sample_col].astype('category').cat.codes + 1
    adata.obs['label_id'] = adata.obs[label_col].map(encode_label_dict).astype(int)
    if batch_col:
        adata.obs['batch_id'] = adata.obs[batch_col].astype('category').cat.codes
    md = MultiAnnData(adata, sampleid='sample_id', dtype=np.float64)
    md.obs_to_sample(["label_id", "batch_id"] if batch_col else ["label_id"])
    # association test
    cna_res = cna.tl.association(md,
                                 getattr(md.samplem, "label_id"),
                                 covs=None,
                                 batches=getattr(md.samplem, "batch_id", None), allow_low_sample_size=True)
    return cna_res, md


