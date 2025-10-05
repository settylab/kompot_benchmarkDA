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
    
