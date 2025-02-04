# cython: language_level=3

import numpy as np

cimport cython
from libc.math cimport sqrt
from libc.math cimport log as std_log
from libc.math cimport exp as std_exp
from libc.math cimport fabs
from libcpp.unordered_map cimport unordered_map

ctypedef unsigned char uint8

ctypedef fused indices_type:
    int
    long

ctypedef fused indptr_type:
    int
    long


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef sgt_cython(const long[:] counts, const long[:] r, const long[:] Nr, const long n0, int size_r, int N, int n_genes):
    cdef Py_ssize_t i
    cdef double x_sum, y_sum, xy_sum, x2_sum, p0, N_star, prob_unseen, sr_cache, sr, sr_next
    cdef long switch_point

    # Calculate Z
    z_buffer = np.zeros(size_r)
    cdef double[:] z = z_buffer

    x_buffer = np.zeros(size_r)
    y_buffer = np.zeros(size_r)
    r_star_buffer = np.zeros(size_r)
    #prob_est_buffer = np.zeros(size_r)
    prop_adj_buffer = np.zeros(n_genes)
    cdef double[:] x = x_buffer
    cdef double[:] y = y_buffer
    cdef double[:] r_star = r_star_buffer
    #cdef double[:] prob_est = prob_est_buffer
    cdef double[:] prop_adj = prop_adj_buffer

    cdef unordered_map[long, double] prob_est

    for i in range(size_r):
        if i == 0:
            z[i] = 2 * Nr[i] / r[i + 1]
        elif i == size_r - 1:
            z[i] = Nr[i] / (r[i] - r[i - 1])
        else:
            z[i] = 2 * Nr[i] / (r[i + 1] - r[i - 1])

    for i in range(size_r):
        x[i] = std_log(r[i])
        y[i] = std_log(z[i])

    # Smooth fit
    x_sum = 0.0
    y_sum = 0.0
    for i in range(size_r):
        x_sum += x[i]
        y_sum += y[i]
        xy_sum += (x[i] * y[i])
        x2_sum += x[i]**2

    slope = (size_r * xy_sum - x_sum * y_sum) / (size_r * x2_sum - x_sum**2)
    intercept = (y_sum - slope * x_sum) / size_r

    # Estimate probability
    p0 = Nr[0] / N if r[0] == 1 else 0
    switch_point = 0
    sr_cache = -1
    for i in range(size_r):
        sr = sr_cache if sr_cache >= 0 else _smoothed_Nr(r[i], slope, intercept)
        sr_next = _smoothed_Nr(r[i] + 1, slope, intercept)
        y_r = (r[i] + 1) * sr_next / sr
        if (i < size_r - 1) and (r[i + 1] == r[i] + 1):
            sr_cache = sr_next    # Cache S(r+1) only when next is r+1
        else:
            sr_cache = -1
        #y_r = (cur_r + 1) * _smoothed_Nr(r[i] + 1, slope, intercept) / _smoothed_Nr(r[i], slope, intercept)
        if switch_point == 0:
            # Decide if switch point
            if (i == size_r - 1) or (r[i + 1] != r[i] + 1):
                switch_point = r[i]
                r_star[i] = y_r
            else:
                x_r = (r[i] + 1) * Nr[i + 1] / Nr[i]
                if fabs(x_r - y_r) <= 1.96 * sqrt((r[i] + 1)**2 * (Nr[i + 1] / Nr[i]**2) * (1 + Nr[i + 1] / Nr[i])):
                    switch_point = r[i]
                    r_star[i] = y_r
                else:
                    r_star[i] = x_r
        else:
            # Compare with switch point
            if r[i] < switch_point:
                r_star[i] = (r[i] + 1) * Nr[i + 1] / Nr[i]
            else:
                r_star[i] = y_r

    N_star = 0.0
    for i in range(size_r):
        N_star += (r_star[i] * Nr[i])
    prob_unseen = p0 / n0 if n0 > 0 else 0
    for i in range(size_r):
        #prob_est[i] = (1 - p0) * r_star[i] / N_star
        prob_est[r[i]] = (1 - p0) * r_star[i] / N_star

    # Assign proportion
    for i in range(n_genes):
        prop_adj[i] = prob_est[counts[i]] if counts[i] > 0 else prob_unseen

    return prop_adj_buffer


cpdef _smoothed_Nr(r, slope, intercept):
    return std_exp(intercept + slope * std_log(r))
