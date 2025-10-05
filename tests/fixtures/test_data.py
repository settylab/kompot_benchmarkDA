"""
Shared test fixtures and sample data for unit and integration tests.
"""

import numpy as np
import pandas as pd
import anndata as ad


def create_mock_adata(n_cells=100, n_genes=50, n_populations=3):
    """
    Create a mock AnnData object for testing.

    Parameters:
    -----------
    n_cells : int
        Number of cells
    n_genes : int
        Number of genes
    n_populations : int
        Number of cell populations

    Returns:
    --------
    adata : AnnData
        Mock AnnData object with populated obs and obsm
    """
    # Create mock expression data
    X = np.random.randn(n_cells, n_genes)

    # Create mock cell metadata
    populations = [f"Pop{(i % n_populations) + 1}" for i in range(n_cells)]
    obs = pd.DataFrame({
        "population": populations,
        "batch": np.random.choice(["B1", "B2"], n_cells),
    }, index=[f"cell_{i}" for i in range(n_cells)])

    # Create mock embeddings
    obsm = {
        "X_pca": np.random.randn(n_cells, 10),
        "DM_EigenVectors": np.random.randn(n_cells, 5),
    }

    adata = ad.AnnData(X=X, obs=obs, obsm=obsm)
    return adata


def create_mock_distances(n_cells=100, n_centroids=3):
    """
    Create mock distance matrix for testing weight calculations.

    Parameters:
    -----------
    n_cells : int
        Number of cells
    n_centroids : int
        Number of centroids

    Returns:
    --------
    distances : ndarray
        Mock distance matrix (n_cells × n_centroids)
    """
    return np.abs(np.random.randn(n_cells, n_centroids))


def create_mock_probabilities(n_cells=100):
    """
    Create mock probability distributions for testing.

    Parameters:
    -----------
    n_cells : int
        Number of cells

    Returns:
    --------
    probs : DataFrame
        Mock probability distributions
    """
    cond1_prob = np.random.uniform(0.3, 0.7, n_cells)
    cond2_prob = 1 - cond1_prob

    return pd.DataFrame({
        "Condition1": cond1_prob,
        "Condition2": cond2_prob
    })
