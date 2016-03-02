import numpy as np

def robustRMS(x):
    """
    Use the inter-quartile range to compute RMS. Effectively rejects outliers.
    """
    iqr = np.percentile(x,75)-np.percentile(x,25)
    rms = iqr/1.349 #approximation
    return rms
