# cython: language_level=3

import numpy as np

cimport cython
from libc.math cimport sqrt


ctypedef unsigned char uint8

ctypedef fused indices_type:
    int
    long

ctypedef fused indptr_type:
    int
    long


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calc_mean_sparse(int M, int N, const float[:] data, indices_type[:] indices, indptr_type[:] indptr, int axis):
    cdef Py_ssize_t i, j

    results = np.zeros(N if axis == 0 else M, dtype = np.float64)
    cdef double[:] res = results

    for i in range(M):
        for j in range(indptr[i], indptr[i + 1]):
            if axis == 0:
                res[indices[j]] += data[j]
            else:
                res[i] += data[j]

    np.divide(results, M if axis == 0 else N, out = results)

    return results


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple calc_mean_and_var_sparse(int M, int N, const float[:] data, indices_type[:] indices, indptr_type[:] indptr, int axis):
    cdef Py_ssize_t i, j

    cdef int size
    cdef double value

    size = N if axis == 0 else M
    mean = np.zeros(size, dtype = np.float64)
    var = np.zeros(size, dtype = np.float64)

    cdef double[:] mean_view = mean
    cdef double[:] var_view = var

    for i in range(M):
        for j in range(indptr[i], indptr[i + 1]):
            value = data[j]
            if axis == 0:
                mean_view[indices[j]] += value
                var_view[indices[j]] += value * value
            else:
                mean_view[i] += value
                var_view[i] += value * value

    size = M if axis == 0 else N
    for i in range(mean_view.size):
        mean_view[i] /= size
        var_view[i] = (var_view[i] - size * mean_view[i] * mean_view[i]) / (size - 1)

    return mean, var


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple calc_mean_and_var_dense(int M, int N, float[:, :] X, int axis):
    cdef Py_ssize_t i, j

    cdef int size
    cdef double value

    size = N if axis == 0 else M
    mean = np.zeros(size, dtype = np.float64)
    var = np.zeros(size, dtype = np.float64)

    cdef double[:] mean_view = mean
    cdef double[:] var_view = var

    for i in range(M):
        for j in range(N):
            value = X[i, j]
            if axis == 0:
                mean_view[j] += value
                var_view[j] += value * value
            else:
                mean_view[i] += value
                var_view[i] += value * value

    size = M if axis == 0 else N
    for i in range(mean_view.size):
        mean_view[i] /= size
        var_view[i] = (var_view[i] - size * mean_view[i] * mean_view[i]) / (size - 1)

    return mean, var


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple calc_stat_per_batch_sparse(int M, int N, const float[:] data, indices_type[:] indices, indptr_type[:] indptr, int nbatch, const int[:] codes):
    cdef Py_ssize_t i, j

    cdef indices_type col
    cdef int code
    cdef double value

    ncells = np.zeros(nbatch, dtype = np.int32)
    means = np.zeros((N, nbatch), dtype = np.float64)
    partial_sum = np.zeros((N, nbatch), dtype = np.float64)

    cdef int[:] ncells_view = ncells
    cdef double[:, :] means_view = means
    cdef double[:, :] ps_view = partial_sum

    for i in range(M):
        code = codes[i]
        ncells_view[code] += 1
        for j in range(indptr[i], indptr[i + 1]):
            col = indices[j]
            value = data[j]
            means_view[col, code] += value
            ps_view[col, code] += value * value

    for j in range(nbatch):
        if ncells_view[j] > 1:
            for i in range(N):
                means_view[i, j] /= ncells_view[j]
                ps_view[i, j] = ps_view[i, j] - ncells_view[j] * means_view[i, j] * means_view[i, j]

    return ncells, means, partial_sum


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple calc_stat_per_batch_dense(int M, int N, const float[:, :] X, int nbatch, const int[:] codes):
    cdef Py_ssize_t i, j

    cdef int code, col
    cdef double value

    ncells = np.zeros(nbatch, dtype = np.int32)
    means = np.zeros((N, nbatch), dtype = np.float64)
    partial_sum = np.zeros((N, nbatch), dtype = np.float64)

    cdef int[:] ncells_view = ncells
    cdef double[:, :] means_view = means
    cdef double[:, :] ps_view = partial_sum

    for i in range(M):
        code = codes[i]
        ncells_view[code] += 1
        for j in range(N):
            value = X[i, j]
            means_view[j, code] += value
            ps_view[j, code] += value * value

    for j in range(nbatch):
        if ncells_view[j] > 1:
            for i in range(N):
                means_view[i, j] /= ncells_view[j]
                ps_view[i, j] = ps_view[i, j] - ncells_view[j] * means_view[i, j] * means_view[i, j]

    return ncells, means, partial_sum


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef normalize_by_count_sparse(int M, int N, float[:] data, indices_type[:] indices, indptr_type[:] indptr, uint8[:] robust, double norm_count):
    cdef Py_ssize_t i, j

    scale = np.zeros(M, dtype = np.float64)
    cdef double[:] scale_factor = scale

    for i in range(M):
        for j in range(indptr[i], indptr[i + 1]):
            if robust[indices[j]]:
                scale_factor[i] += data[j]
        scale_factor[i] = norm_count / scale_factor[i]
        for j in range(indptr[i], indptr[i + 1]):
            data[j] *= scale_factor[i]

    return scale


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef normalize_by_count_dense(int M, int N, float[:, :] X, uint8[:] robust, double norm_count):
    cdef Py_ssize_t i, j

    scale = np.zeros(M, dtype = np.float64)
    cdef double[:] scale_factor = scale

    for i in range(M):
        for j in range(N):
            if robust[j]:
                scale_factor[i] += X[i, j]
        scale_factor[i] = norm_count / scale_factor[i]
        for j in range(N):
            X[i, j] *= scale_factor[i]

    return scale


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calc_sig_background_sparse(int M, int N, const float[:] data, indices_type[:] indices, indptr_type[:] indptr, int n_bins, const int[:] codes, const double[:] mean_vec):
    cdef Py_ssize_t i, j

    sig_bkg_mean = np.zeros((M, n_bins), dtype = np.float64)
    sig_bkg_std = np.ones((M, n_bins), dtype = np.float64)
    bin_mean = np.zeros(n_bins, dtype = np.float64)
    bin_sq = np.zeros(n_bins, dtype = np.float64)
    cdef double[:, :] sig_bkgm_view = sig_bkg_mean
    cdef double[:, :] sig_bkgs_view = sig_bkg_std
    cdef double[:] bin_mean_view = bin_mean
    cdef double[:] bin_sq_view = bin_sq

    cdef int[:] bin_count = np.zeros(n_bins, dtype = np.int32)
    cdef double estd

    for j in range(N):
        bin_count[codes[j]] += 1
        bin_mean_view[codes[j]] += mean_vec[j]
        bin_sq_view[codes[j]] += (mean_vec[j] * mean_vec[j])

    for i in range(M):
        for j in range(indptr[i], indptr[i + 1]):
            sig_bkgm_view[i, codes[indices[j]]] += data[j]
            sig_bkgs_view[i, codes[indices[j]]] += (data[j]**2 - 2 * data[j] * mean_vec[indices[j]])

    for j in range(n_bins):
        for i in range(M):
            sig_bkgm_view[i, j] = (sig_bkgm_view[i, j] - bin_mean_view[j]) / bin_count[j]
            sig_bkgs_view[i, j] += bin_sq_view[j]
            estd = 0.0
            if bin_count[j] > 1:
                estd = sqrt((sig_bkgs_view[i, j] - bin_count[j] * sig_bkgm_view[i, j] * sig_bkgm_view[i, j]) / (bin_count[j] - 1))
            sig_bkgs_view[i, j] = estd if estd > 1e-4 else 1.0

    return sig_bkg_mean, sig_bkg_std


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calc_sig_background_sparse_no_std(int M, int N, const float[:] data, indices_type[:] indices, indptr_type[:] indptr, int n_bins, const int[:] codes, const double[:] mean_vec):
    cdef Py_ssize_t i, j

    sig_bkg_mean = np.zeros((M, n_bins), dtype = np.float64)
    bin_mean = np.zeros(n_bins, dtype = np.float64)
    cdef double[:, :] sig_bkgm_view = sig_bkg_mean
    cdef double[:] bin_mean_view = bin_mean

    cdef int[:] bin_count = np.zeros(n_bins, dtype = np.int32)

    for j in range(N):
        bin_count[codes[j]] += 1
        bin_mean_view[codes[j]] += mean_vec[j]

    for i in range(M):
        for j in range(indptr[i], indptr[i + 1]):
            sig_bkgm_view[i, codes[indices[j]]] += data[j]

    for j in range(n_bins):
        for i in range(M):
            sig_bkgm_view[i, j] = (sig_bkgm_view[i, j] - bin_mean_view[j]) / bin_count[j]

    return sig_bkg_mean


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calc_sig_background_dense(int M, int N, const float[:, :] X, int n_bins, const int[:] codes, const double[:] mean_vec):
    cdef Py_ssize_t i, j

    sig_bkg_mean = np.zeros((M, n_bins), dtype = np.float64)
    sig_bkg_std = np.ones((M, n_bins), dtype = np.float64)
    cdef double[:, :] sig_bkgm_view = sig_bkg_mean
    cdef double[:, :] sig_bkgs_view = sig_bkg_std

    cdef int[:] bin_count = np.zeros(n_bins, dtype = np.int32)
    cdef double cexpr, estd

    for j in range(N):
        bin_count[codes[j]] += 1

    for i in range(M):
        for j in range(N):
            cexpr = X[i, j] - mean_vec[j]
            sig_bkgm_view[i, codes[j]] += cexpr
            sig_bkgs_view[i, codes[j]] += cexpr * cexpr

    for j in range(n_bins):
        for i in range(M):
            sig_bkgm_view[i, j] /= bin_count[j]
            estd = 0.0
            if bin_count[j] > 1:
                estd = sqrt((sig_bkgs_view[i, j] - bin_count[j] * sig_bkgm_view[i, j] * sig_bkgm_view[i, j]) / (bin_count[j] - 1))
            sig_bkgs_view[i, j] = estd if estd > 1e-4 else 1.0

    return sig_bkg_mean, sig_bkg_std


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calc_sig_background_dense_no_std(int M, int N, const float[:, :] X, int n_bins, const int[:] codes, const double[:] mean_vec):
    cdef Py_ssize_t i, j

    sig_bkg_mean = np.zeros((M, n_bins), dtype = np.float64)
    cdef double[:, :] sig_bkgm_view = sig_bkg_mean

    cdef int[:] bin_count = np.zeros(n_bins, dtype = np.int32)
    cdef double cexpr

    for j in range(N):
        bin_count[codes[j]] += 1

    for i in range(M):
        for j in range(N):
            cexpr = X[i, j] - mean_vec[j]
            sig_bkgm_view[i, codes[j]] += cexpr

    for j in range(n_bins):
        for i in range(M):
            sig_bkgm_view[i, j] /= bin_count[j]

    return sig_bkg_mean


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple simulate_doublets_sparse(int n_sim, int N, const int[:] data, indices_type[:] indices, indptr_type[:] indptr, const int[:, :] doublet_indices):
    ### Assume X is ordered in ascending order
    cdef Py_ssize_t i, j, k, u, v, u_up, v_up, size, counter

    # Generate new count matrix
    size = 0
    for i in range(n_sim):
        u = doublet_indices[i, 0]
        v = doublet_indices[i, 1]
        size += (indptr[u + 1] - indptr[u]) + (indptr[v + 1] - indptr[v])

    indices_dtype = np.int64 if indices_type is long else np.int32
    indptr_dtype = np.int64 if indptr_type is long else np.int32

    out_data_buffer = np.zeros(size, dtype = np.int32)
    out_indices_buffer = np.zeros(size, dtype = indices_dtype)
    out_indptr_buffer = np.zeros(n_sim + 1, dtype = indptr_dtype)

    cdef int[:] out_data = out_data_buffer
    cdef indices_type[:] out_indices = out_indices_buffer
    cdef indptr_type[:] out_indptr = out_indptr_buffer

    counter = 0
    for i in range(n_sim):
        out_indptr[i] = counter

        u = doublet_indices[i, 0]
        v = doublet_indices[i, 1]

        j = indptr[u]
        u_up = indptr[u + 1]
        k = indptr[v]
        v_up = indptr[v + 1]

        while j < u_up and k < v_up:
            if indices[j] < indices[k]:
                out_indices[counter] = indices[j]
                out_data[counter] = data[j]
                j += 1
            elif indices[j] == indices[k]:
                out_indices[counter] = indices[j]
                out_data[counter] = data[j] + data[k]
                j += 1
                k += 1
            else:
                out_indices[counter] = indices[k]
                out_data[counter] = data[k]
                k += 1
            counter += 1

        while j < u_up:
            out_indices[counter] = indices[j]
            out_data[counter] = data[j]
            counter += 1
            j += 1

        while k < v_up:
            out_indices[counter] = indices[k]
            out_data[counter] = data[k]
            counter += 1
            k += 1

    out_indptr[n_sim] = counter
    out_data_buffer.resize(out_indptr[n_sim], refcheck=False)
    out_indices_buffer.resize(out_indptr[n_sim], refcheck=False)

    return out_data_buffer, out_indices_buffer, out_indptr_buffer


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef simulate_doublets_dense(int n_sim, int N, const int[:, :] X, const int[:, :] doublet_indices):
    cdef Py_ssize_t i, j, u, v

    # Generate new count matrix
    results = np.zeros((n_sim, N), dtype = np.int32)
    cdef int[:, :] out_array = results

    for i in range(n_sim):
        u = doublet_indices[i, 0]
        v = doublet_indices[i, 1]
        for j in range(N):
            out_array[i, j] = X[u, j] + X[v, j]

    return results
