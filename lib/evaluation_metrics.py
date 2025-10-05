"""
Centralized evaluation metrics for differential abundance analysis.
All methods should use these standardized implementations.
"""

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score


def modified_auroc(adata, lfc):
    """
    Calculate modified AUROC for DA results.

    Computes separate AUROC scores for positive and negative log fold changes,
    then returns their mean. Does NOT flip scores - poor predictions should
    remain poor (negative scores are valid).

    Parameters:
    -----------
    adata : AnnData
        AnnData object with 'true_labels' in obs
    lfc : array-like
        Log fold change or likelihood scores

    Returns:
    --------
    tuple : (auc_score, auc_score_neg, auc_score_pos)
        Mean AUROC and separate scores for negative and positive classes
    """
    true_label = adata.obs["true_labels"]
    true_label_transformed_neg = (true_label == "NegLFC")
    true_label_transformed_pos = (true_label == "PosLFC")

    # Compute AUROC for negative class
    if np.any(true_label_transformed_neg) and np.any(~true_label_transformed_neg):
        auc_score_neg = roc_auc_score(true_label_transformed_neg, -lfc)
    else:
        auc_score_neg = np.nan

    # Compute AUROC for positive class
    if np.any(true_label_transformed_pos) and np.any(~true_label_transformed_pos):
        auc_score_pos = roc_auc_score(true_label_transformed_pos, lfc)
    else:
        auc_score_pos = np.nan

    auc_score = np.nanmean([auc_score_neg, auc_score_pos])
    return auc_score, auc_score_neg, auc_score_pos


def modified_auprc(adata, lfc):
    """
    Calculate modified AUPRC for DA results.

    Computes separate AUPRC scores for positive and negative log fold changes,
    then returns their mean.

    Parameters:
    -----------
    adata : AnnData
        AnnData object with 'true_labels' in obs
    lfc : array-like
        Log fold change or likelihood scores

    Returns:
    --------
    tuple : (auprc_score, auprc_score_neg, auprc_score_pos)
        Mean AUPRC and separate scores for negative and positive classes
    """
    true_label = adata.obs["true_labels"]
    true_label_transformed_neg = (true_label == "NegLFC")
    true_label_transformed_pos = (true_label == "PosLFC")

    # Compute AUPRC for negative class
    if np.any(true_label_transformed_neg) and np.any(~true_label_transformed_neg):
        auprc_score_neg = average_precision_score(true_label_transformed_neg, -lfc)
    else:
        auprc_score_neg = np.nan

    # Compute AUPRC for positive class
    if np.any(true_label_transformed_pos) and np.any(~true_label_transformed_pos):
        auprc_score_pos = average_precision_score(true_label_transformed_pos, lfc)
    else:
        auprc_score_pos = np.nan

    auprc_score = np.nanmean([auprc_score_neg, auprc_score_pos])
    return auprc_score, auprc_score_neg, auprc_score_pos


def calculate_outcome(adata, true_label, predicted_label):
    """
    Calculate classification metrics (TP, FP, TN, FN, etc.) for DA predictions.

    Parameters:
    -----------
    adata : AnnData
        AnnData object with obs containing the label columns
    true_label : str
        Column name in adata.obs for true labels
    predicted_label : str
        Column name in adata.obs for predicted labels

    Returns:
    --------
    pd.DataFrame
        DataFrame with metrics: TP, FP, FN, TN, TPR, FPR, TNR, FNR,
        FDR, Precision, Power, Accuracy
    """
    TP = ((adata.obs[true_label] == adata.obs[predicted_label]) &
          (adata.obs[predicted_label] != 'NotDA')).sum()
    FP = ((adata.obs[true_label] != adata.obs[predicted_label]) &
          (adata.obs[predicted_label] != 'NotDA')).sum()
    FN = ((adata.obs[true_label] != adata.obs[predicted_label]) &
          (adata.obs[predicted_label] == 'NotDA')).sum()
    TN = ((adata.obs[true_label] == adata.obs[predicted_label]) &
          (adata.obs[predicted_label] == 'NotDA')).sum()

    metrics_dic = {
        'TP': TP, 'FP': FP,
        'FN': FN, 'TN': TN,
        'TPR': [TP / (TP + FN)] if (TP + FN) > 0 else [np.nan],
        'FPR': [FP / (FP + TN)] if (FP + TN) > 0 else [np.nan],
        'TNR': [TN / (TN + FP)] if (TN + FP) > 0 else [np.nan],
        'FNR': [FN / (FN + TP)] if (FN + TP) > 0 else [np.nan],
        'FDR': [FP / (TP + FP)] if (TP + FP) > 0 else [np.nan],
        'Precision': [TP / (TP + FP)] if (TP + FP) > 0 else [np.nan],
        'Power': [1 - FN / (FN + TP)] if (FN + TP) > 0 else [np.nan],
        'Accuracy': [(TP + TN) / (TP + TN + FP + FN)] if (TP + TN + FP + FN) > 0 else [np.nan]
    }

    return pd.DataFrame(metrics_dic, index=['metric'])
