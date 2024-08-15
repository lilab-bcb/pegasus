import numpy as np
from scipy.stats import linregress


def simple_good_turing(counts):
    r, Nr = np.unique(counts, return_counts=True)
    if r[0] == 0:
        N0 = Nr[0]
        Nr = Nr[1:]
        r = r[1:]
    else:
        N0 = 0

    # Perform linear regression on log-log scale
    log_r = np.log(r)
    log_Nr = np.log(Nr)
    slope, intercept, _, _, _ = linregress(log_r, log_Nr)

    # Calculate smoothed frequencies
    smoothed_Nr = np.exp(intercept + slope * log_r)

    # Calculate Good-Turing estimates
    gt_estimates = {}
    for i, freq in enumerate(r):
        if i < len(r) - 1:
            gt_estimates[freq] = ((freq + 1) * np.exp(intercept + slope * np.log(freq + 1))) / np.exp(intercept + slope * np.log(freq))
        else:
            gt_estimates[freq] = freq

    # Calculate the total number of items
    #N = np.sum(r * Nr)
    N = np.sum(r * smoothed_Nr)

    # Calculate the probability for unseen items
    if N0 != 0:
        p0 = N0 / N
    elif 1 in r:
        p0 =  Nr[np.where(r==1)[0]] / N
    else:
        p0 = 0

    # Calculate the adjusted probabilities for seen items
    adjusted_probs = {}
    for freq in r:
        adjusted_probs[freq] = gt_estimates[freq] / N

    return adjusted_probs, p0
