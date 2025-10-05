import numpy as np
from sklearn.preprocessing import StandardScaler
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.metrics import precision_recall_curve, average_precision_score


def scale(x, min_x, max_x):
    """
    Linearly rescale array values to a new range.

    Maps the range [min(x), max(x)] to [min_x, max_x] using linear interpolation.
    Preserves relative differences between values while changing the scale.

    Parameters:
    - x: Array-like of values to rescale
    - min_x: Target minimum value
    - max_x: Target maximum value

    Returns:
    - scaled_x: Array with values rescaled to [min_x, max_x]

    Example:
        scale([0, 0.5, 1.0], 0.5, 0.95)
        → [0.5, 0.725, 0.95]
    """
    scaled_x = (((x - np.min(x)) / (np.max(x) - np.min(x))) * (max_x - min_x)) + min_x
    return scaled_x


def log_function(x, a):
    """
    Apply sigmoid (logistic) function to transform weight matrix.

    Sigmoid function: f(x) = 1 / (1 + exp(-a*x))
    Maps values to (0, 1) range with smooth S-shaped curve.

    Parameters:
    - x: Array-like of input values (typically standardized weights)
    - a: Slope parameter controlling steepness of sigmoid
         (larger a → steeper curve, smaller a → gentler curve)

    Returns:
    - Sigmoid-transformed values in (0, 1) range

    Note: Input x should be standardized (mean=0, std=1) for consistent behavior.
    """
    log_data_denominator = 1 + np.exp(-a * x)
    log_data = 1 / (log_data_denominator)
    return log_data


def normalization(w):
    """
    Standardize weight matrix using z-score normalization.

    Transforms each column to have mean=0 and std=1 using sklearn's StandardScaler.
    This ensures consistent scale across different distance metrics.

    Parameters:
    - w: Array-like weight matrix (cells x populations)

    Returns:
    - w_scaled: Standardized weight matrix (mean=0, std=1 per column)
    """
    scaler = StandardScaler()
    w_scaled = scaler.fit_transform(w)
    return w_scaled


def convert_number_str(number_str):
    """
    Convert numeric string to int or float depending on value.

    If the number is a whole number (e.g., "1.0"), returns int.
    Otherwise returns float.

    Parameters:
    - number_str: String representation of a number

    Returns:
    - int if whole number, float otherwise

    Example:
        convert_number_str("1.0") → 1
        convert_number_str("1.5") → 1.5
    """
    number_float = float(number_str)

    if number_float.is_integer():
        return int(number_float)
    else:
        return number_float
