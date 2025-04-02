import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score

def scale(x,min_x, max_x):
    
    """
    scale condition probability
    """
    scaled_x = (((x-np.min(x))/(np.max(x)-np.min(x)))*(max_x - min_x)) + min_x
    return scaled_x

def log_function(x,a):
    """
    sigmoid function for normalized weight matrix
    """
    log_data_denominator = 1+np.exp(-a * x)
    log_data = 1/(log_data_denominator)
    return log_data

def normalization(w):
    """
    Used sklearn function StandardScaler to normalize weight matrix
    """
    scaler = StandardScaler()
    w_scaled = scaler.fit_transform(w)
    return w_scaled


def convert_number_str(number_str):
    # Step 1: Convert the string to a float
    number_float = float(number_str)
    
    # Step 2: Check if the float is a whole number
    if number_float.is_integer():
        # Step 3: Convert to an integer if it is a whole number
        return int(number_float)
    else:
        # Keep it as a float if it is not a whole number
        return number_float
    


def modified_auroc(adata, lfc):
    true_label = adata.obs["true_labels"]
    true_label_transformed_neg = (true_label == "NegLFC").astype(int)
    true_label_transformed_pos = (true_label == "PosLFC").astype(int)
    
    try:
        if np.any(true_label_transformed_neg) and np.any(~true_label_transformed_neg):
            auc_score_neg = roc_auc_score(true_label_transformed_neg, lfc)
            if auc_score_neg <= 0.5:
                auc_score_neg = 1 - auc_score_neg
        else:
            auc_score_neg = np.nan
        
        if np.any(true_label_transformed_pos) and np.any(~true_label_transformed_pos):
            auc_score_pos = roc_auc_score(true_label_transformed_pos, lfc)
            if auc_score_pos <= 0.5:
                auc_score_pos = 1 - auc_score_pos
        else:
            auc_score_pos = np.nan
    except ValueError as e:
        print(f"ValueError: {e}")
        auc_score_neg = np.nan
        auc_score_pos = np.nan
    
    auc_score = np.nanmean([auc_score_neg, auc_score_pos])
    return auc_score

def modified_auprc(adata, lfc):
    true_label = adata.obs["true_labels"]
    true_label_transformed_neg = (true_label == "NegLFC").astype(int)
    true_label_transformed_pos = (true_label == "PosLFC").astype(int)
    
    try:
        # Check if there are both positive and negative samples in the transformed arrays
        if np.any(true_label_transformed_neg) and np.any(~true_label_transformed_neg):
            auprc_score_neg = average_precision_score(true_label_transformed_neg, lfc)
        else:
            auprc_score_neg = np.nan
        
        if np.any(true_label_transformed_pos) and np.any(~true_label_transformed_pos):
            auprc_score_pos = average_precision_score(true_label_transformed_pos, lfc)
        else:
            auprc_score_pos = np.nan
    except ValueError as e:
        print(f"ValueError: {e}")
        auprc_score_neg = np.nan
        auprc_score_pos = np.nan
    
    auprc_score = np.nanmean([auprc_score_neg, auprc_score_pos])
    return auprc_score