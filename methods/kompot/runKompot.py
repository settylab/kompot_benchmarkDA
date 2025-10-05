import kompot
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import palantir
import logging
import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from lib import shared_embedding_utils

logger = logging.getLogger("kompot")

def runKOMPOT(
    adata,
    label_col: str,
    obsm_key: str = "X_pca",
    ls_factor: float = 10.0,
    n_landmarks: int = None,
    random_state: int = None,
    **kwargs
):
    """
    Run Kompot differential abundance analysis.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object containing the dataset
    label_col : str
        Column name in adata.obs containing condition labels
    obsm_key : str, default "X_pca"
        Key in adata.obsm to use for embedding (e.g., "X_pca", "DM_EigenVectors")
        The embedding MUST already exist in adata.obsm - this function will not compute it.
    ls_factor : float, default 10.0
        Length scale factor for density estimation
    n_landmarks : int, optional
        Number of landmarks to use. If None, kompot will determine automatically
    random_state : int, optional
        Random seed for reproducibility
    **kwargs : additional arguments
        Additional keyword arguments passed to kompot.compute_differential_abundance

    Returns:
    --------
    tuple
        (log_fold_change, z_scores) arrays
    """

    # Get unique conditions
    conditions = adata.obs[label_col].unique()
    if len(conditions) != 2:
        raise ValueError(f"Expected exactly 2 conditions, got {len(conditions)}: {conditions}")

    condition1, condition2 = sorted(conditions)

    logger.info(f"Running Kompot DA analysis comparing {condition1} vs {condition2}")
    logger.info(f"Using embedding: {obsm_key}")
    logger.info(f"Number of cells: {adata.n_obs}")

    # Verify the obsm_key exists - fail if it doesn't
    if obsm_key not in adata.obsm:
        raise ValueError(
            f"Embedding '{obsm_key}' not found in adata.obsm. "
            f"Available embeddings: {list(adata.obsm.keys())}. "
            f"Use runKOMPOT_with_params() or ensure embeddings are computed before calling this function."
        )

    try:
        # Run kompot differential abundance analysis
        results = kompot.compute_differential_abundance(
            adata,
            groupby=label_col,
            condition1=condition1,
            condition2=condition2,
            obsm_key=obsm_key,
            ls_factor=ls_factor,
            n_landmarks=n_landmarks,
            random_state=random_state,
            return_full_results=True,
            inplace=False,
            **kwargs
        )
        log_fold_change = np.asarray(results["log_fold_change"])
        z_scores = np.asarray(results['log_fold_change_zscore'])

        logger.info(f"Kompot analysis completed successfully")
        logger.info(f"Log fold change range: [{np.min(log_fold_change):.3f}, {np.max(log_fold_change):.3f}]")
        logger.info(f"Z-scores range: [{np.min(z_scores):.3f}, {np.max(z_scores):.3f}]")

        return log_fold_change, z_scores

    except Exception as e:
        logger.error(f"Error running Kompot: {str(e)}")
        raise


def runKOMPOT_with_params(
    adata,
    label_col: str,
    use_dm: bool = True,
    dm_comp: int = 10,
    ls_factor: float = 10.0,
    **kwargs
):
    """
    Wrapper function that ensures consistent embeddings across all methods.

    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    label_col : str
        Column name for condition labels
    use_dm : bool, default True
        Whether to use diffusion maps (True) or batch-corrected PCA (False)
    dm_comp : int, default 10
        Number of diffusion map components
    ls_factor : float, default 10.0
        Kompot's default length scale factor
    **kwargs : additional arguments
        Passed to runKOMPOT
    """

    # Ensure consistent embeddings using shared utility
    adata = shared_embedding_utils.ensure_batch_corrected_embeddings(
        adata, dm_comp=dm_comp
    )

    # Get the appropriate obsm key
    obsm_key = shared_embedding_utils.get_obsm_key_for_method(use_dm=use_dm)

    logger.info(f"Kompot using embedding: {obsm_key}")
    if use_dm:
        logger.info(f"Diffusion maps with {dm_comp} components (computed from batch-corrected PCA)")
    else:
        logger.info("Batch-corrected PCA")

    return runKOMPOT(
        adata=adata,
        label_col=label_col,
        obsm_key=obsm_key,
        ls_factor=ls_factor,
        **kwargs
    )


def threshold_kompot(neg_log10_ptp_scores):
    """
    Generate threshold range for kompot neg_log10_ptp scores.
    Similar to threshold_mellon but for kompot's significance scores.

    Parameters:
    -----------
    neg_log10_ptp_scores : np.ndarray
        Negative log10 PTP scores from kompot

    Returns:
    --------
    np.ndarray
        Array of thresholds for classification
    """
    lower = 0
    upper = np.max(np.abs(neg_log10_ptp_scores)) - 1e-8
    assert lower < upper, "lower bound is not less than upper bound"
    thresholds = np.arange(lower, upper, (upper - lower) / 100)
    return thresholds


def kompot2output(lfcs, significance_scores, out_type="continuous", thresholds=None):
    """
    Convert kompot results to output format for benchmarking.
    Similar to mellon2output but for kompot results.

    Parameters:
    -----------
    lfcs : np.ndarray
        Log fold change values
    significance_scores : np.ndarray
        Signed significance scores (neg_log10_ptp * sign(lfc))
    out_type : str, default "continuous"
        Output type, either "continuous" or classification
    thresholds : float, list, or np.ndarray, optional
        Threshold(s) for classification

    Returns:
    --------
    np.ndarray or list
        If out_type="continuous", returns significance_scores
        Otherwise returns classified labels
    """
    if out_type == "continuous":
        da_cell = significance_scores
    else:
        def get_da_cell(thres):
            issign = np.abs(significance_scores) > thres
            isPos = np.logical_and(issign, lfcs > 0)
            isNeg = np.logical_and(issign, lfcs < 0)

            da = np.array(["NotDA"] * len(significance_scores), dtype=object)
            da[isPos] = "PosLFC"
            da[isNeg] = "NegLFC"
            return da

        if isinstance(thresholds, float):
            da_cell = get_da_cell(thresholds)
        elif isinstance(thresholds, list) or isinstance(thresholds, np.ndarray):
            da_cell = []
            for _, thres in enumerate(thresholds):
                da_cell.append(get_da_cell(thres))
        else:
            raise RuntimeError("param: thresholds can only support list or float")

    return da_cell


