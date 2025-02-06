import numpy as np

from numba import njit
from numba.core import types
from numba.typed import Dict

from pegasus.cylib.fast_sgt import sgt_cython

def sgt_estimate_cython(counts):
    r, Nr = np.unique(counts, return_counts=True)
    if r[0] == 0:
        n0 = Nr[0]
        Nr = Nr[1:].view()
        r = r[1:].view()
    else:
        n0 = 0

    N = np.sum(counts)
    prop_adj = sgt_cython(counts, r, Nr, n0, r.size, N, counts.size)
    return prop_adj


def sgt_estimate(
    counts: np.array,
) -> np.array:
    """Perform Simple Good-Turing (SGT) estimation on proportions of each species/event from an observed sample.
    The estimation uses Simple Good-Turing method ([Gale95]_), which has the nice property of assigning non-zero proportions
    to species/events not observed in the given sample.
    An example use case is to estimate proportion of each gene in ambient, given UMI counts of each gene from a selected ambient sample.

    Parameters
    -----------
    counts: ``numpy.array``
        The 1D array of shape ``(n_species,)``, where each element is an integer count of the corresponding species/event.

    Returns
    -------
    An array of shape ``(n_species,)``, where each element is a proportion estimated by SGT for the species/event.

    Examples
    ---------
    >>> pg.sgt_estimate(counts)
    """
    r, Nr = np.unique(counts, return_counts=True)
    if r[0] == 0:
        n0 = Nr[0]
        Nr = Nr[1:].view()
        r = r[1:].view()
    else:
        n0 = 0

    # Smooth fit by simple linear regression
    z = _calc_z(r, Nr)
    x = np.log(r)
    y = np.log(z)
    slope, intercept = _smooth_fit(x, y)

    prop_adj = _estimate_proportion(counts, r, Nr, n0, slope, intercept)
    return prop_adj


@njit
def _estimate_proportion(counts, r, Nr, n0, slope, intercept):
    N = np.sum(r * Nr)
    p0 = Nr[0] / N if r[0] == 1 else 0
    switch_point = 0
    r_star = np.zeros(r.size)

    sr_cache = np.nan
    for i, cur_r in enumerate(r):
        sr = sr_cache if not np.isnan(sr_cache) else _smoothed_Nr(cur_r, slope, intercept)
        sr_next = _smoothed_Nr(cur_r + 1, slope, intercept)
        y_r = (cur_r + 1) * sr_next / sr
        if (i < r.size - 1) and (r[i + 1] == cur_r + 1):
            sr_cache = sr_next  # Cache S(r+1) only when next is r+1
        else:
            sr_cache = np.nan
        if switch_point == 0:
            # Decide if switch point
            if (i == r.size - 1) or (r[i + 1] != cur_r + 1):
                switch_point = cur_r
                r_star[i] = y_r
            else:
                x_r = (cur_r + 1) * Nr[i + 1] / Nr[i]
                if np.abs(x_r - y_r) <= 1.96 * np.sqrt((cur_r + 1)**2 * (Nr[i + 1] / Nr[i]**2) * (1 + Nr[i + 1] / Nr[i])):
                    switch_point = cur_r
                    r_star[i] = y_r
                else:
                    r_star[i] = x_r
        else:
            # Compare with switch point
            if cur_r < switch_point:
                r_star[i] = (cur_r + 1) * Nr[i + 1] / Nr[i]
            else:
                r_star[i] = y_r

    N_star = np.sum(r_star * Nr)
    prob_unseen = p0 / n0 if n0 > 0 else 0
    prob_est = Dict.empty(
        key_type=types.int64,
        value_type=types.float64,
    )
    for i, cur_r in enumerate(r):
        prob_est[cur_r] = (1 - p0) * r_star[i] / N_star

    prop_adj = np.full(counts.size, np.nan)
    for i in np.arange(counts.size):
        prop_adj[i] = prob_est[counts[i]] if counts[i] > 0 else prob_unseen

    return prop_adj


@njit
def _calc_z(r, Nr):
    z = np.zeros(r.size)
    for i in np.arange(z.size):
        if i == 0:
            z[i] = 2 * Nr[i] / r[i + 1]
        elif i == z.size - 1:
            z[i] = Nr[i] / (r[i] - r[i - 1])
        else:
            z[i] = 2 * Nr[i] / (r[i + 1] - r[i - 1])
    return z


@njit
def _smooth_fit(x, y):
    n_obs = x.size
    x_sum = x.sum()
    y_sum = y.sum()

    slope = (n_obs * np.sum(x*y) - x_sum * y_sum) / (n_obs * np.sum(x**2) - x_sum**2)
    intercept = (y_sum - slope * x_sum) / n_obs

    return slope, intercept


@njit
def _smoothed_Nr(r, slope, intercept):
    return np.exp(intercept + slope * np.log(r))
