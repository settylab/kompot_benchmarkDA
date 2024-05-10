import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd


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