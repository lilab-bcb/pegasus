# cython: language_level=3

import numpy as np
from scipy.sparse import csr_matrix, issparse
from typing import Union, List

cimport cython

ctypedef unsigned char uint8



@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calc_mean(X: Union[csr_matrix, np.ndarray], axis: int):
    if not issparse(X):
        return X.mean(axis = axis, dtype = np.float64)

    cdef Py_ssize_t i, j
    cdef int M = X.shape[0] # M rows
    cdef int N = X.shape[1] # N columns

    cdef const float[:] data = X.data
    cdef const int[:] indices = X.indices
    cdef const long[:] indptr = X.indptr.astype(np.int64)

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
cpdef tuple calc_mean_and_var(X: Union[csr_matrix, np.ndarray], axis: int):
    cdef Py_ssize_t i, j

    cdef int M = X.shape[0] # M rows
    cdef int N = X.shape[1] # N columns
    cdef int size

    # if sparse
    cdef const float[:] data
    cdef const int[:] indices
    cdef const long[:] indptr

    # if not sparse
    cdef const float[:, :] data_array 

    
    cdef double value

    size = N if axis == 0 else M
    mean = np.zeros(size, dtype = np.float64)
    var = np.zeros(size, dtype = np.float64)

    cdef double[:] mean_view = mean
    cdef double[:] var_view = var

    if issparse(X):
        data = X.data
        indices = X.indices
        indptr = X.indptr.astype(np.int64)

        for i in range(M):
            for j in range(indptr[i], indptr[i + 1]):
                value = data[j]
                if axis == 0:
                    mean_view[indices[j]] += value
                    var_view[indices[j]] += value * value
                else:
                    mean_view[i] += value
                    var_view[i] += value * value
    else:
        data_array = X

        for i in range(M):
            for j in range(N):
                value = data_array[i, j]
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
cpdef tuple calc_stat_per_batch(X: Union[csr_matrix, np.ndarray], batch):
    cdef Py_ssize_t i, j

    cdef int M = X.shape[0] # M rows
    cdef int N = X.shape[1] # N columns

    # if sparse
    cdef const float[:] data
    cdef const int[:] indices
    cdef const long[:] indptr

    # if not sparse
    cdef const float[:, :] data_array 

    # for batch
    cdef int nbatch = batch.categories.size
    cdef const int[:] codes = batch.codes.astype(np.int32)

    cdef int code, col
    cdef double value

    ncells = np.zeros(nbatch, dtype = np.int32)
    means = np.zeros((N, nbatch), dtype = np.float64)
    partial_sum = np.zeros((N, nbatch), dtype = np.float64)

    cdef int[:] ncells_view = ncells
    cdef double[:, :] means_view = means
    cdef double[:, :] ps_view = partial_sum

    if issparse(X):
        data = X.data
        indices = X.indices
        indptr = X.indptr.astype(np.int64)

        for i in range(M):
            code = codes[i]
            ncells_view[code] += 1
            for j in range(indptr[i], indptr[i + 1]):
                col = indices[j]
                value = data[j]
                means_view[col, code] += value
                ps_view[col, code] += value * value
    else:
        data_array = X

        for i in range(M):
            code = codes[i]
            ncells_view[code] += 1
            for j in range(N):
                value = data_array[i, j]
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
cpdef normalize_by_count(X: Union[csr_matrix, np.ndarray], uint8[:] robust, double norm_count, bint log_transform):
    cdef Py_ssize_t i, j

    cdef int M = X.shape[0] # M rows
    cdef int N = X.shape[1] # N columns

    # if sparse
    cdef float[:] data
    cdef const int[:] indices
    cdef const long[:] indptr

    # if not sparse
    cdef float[:, :] data_array 

    scale = np.zeros(M, dtype = np.float64)
    cdef double[:] scale_factor = scale

    if issparse(X):
        data = X.data
        indices = X.indices
        indptr = X.indptr.astype(np.int64)

        for i in range(M):
            for j in range(indptr[i], indptr[i + 1]):
                if robust[indices[j]]:
                    scale_factor[i] += data[j]
            scale_factor[i] = norm_count / scale_factor[i]
            for j in range(indptr[i], indptr[i + 1]):
                data[j] *= scale_factor[i]

        if log_transform:
            np.log1p(X.data, out = X.data)
    else:
        data_array = X

        for i in range(M):
            for j in range(N):
                if robust[j]:
                    scale_factor[i] += data_array[i, j]
            scale_factor[i] = norm_count / scale_factor[i]
            for j in range(N):
                data_array[i, j] *= scale_factor[i]

        if log_transform:
            np.log1p(X, out = X)

    return scale


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef calc_sig_background(X: Union[csr_matrix, np.ndarray], bins, const double[:] mean_vec):
    cdef Py_ssize_t i, j

    cdef int M = X.shape[0] # M rows
    cdef int N = X.shape[1] # N columns

    # if sparse
    cdef const float[:] data
    cdef const int[:] indices
    cdef const long[:] indptr

    # if not sparse
    cdef const float[:, :] data_array

    # for bins
    cdef int n_bins = bins.categories.size
    cdef const int[:] codes = bins.codes.astype(np.int32)

    sig_background = np.zeros((M, n_bins), dtype = np.float64)
    cdef double[:, :] sig_back_view = sig_background

    cdef int[:] bin_count = np.zeros(n_bins, dtype = np.int32)
    cdef double[:] bin_avg = np.zeros(n_bins, dtype = np.float64)

    for j in range(N):
        bin_count[codes[j]] += 1
        bin_avg[codes[j]] += mean_vec[j]

    if issparse(X):
        data = X.data
        indices = X.indices
        indptr = X.indptr.astype(np.int64)

        for i in range(M):
            for j in range(indptr[i], indptr[i + 1]):
                sig_back_view[i, codes[indices[j]]] += data[j]
    else:
        data_array = X
        for i in range(M):
            for j in range(N):
                sig_back_view[i, codes[j]] += data_array[i, j]

    for j in range(n_bins):
        bin_avg[j] /= bin_count[j]
        for i in range(M):
            sig_back_view[i, j] = sig_back_view[i, j] / bin_count[j] - bin_avg[j]

    return sig_background


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple simulate_doublets(X: Union[csr_matrix, np.ndarray], double sim_doublet_ratio, Py_ssize_t random_state = 0):
    ### Also return simulated index in case there is a need
    ### Assume X is ordered in ascending order
    cdef Py_ssize_t i, j, k, u, v, u_up, v_up, size, counter

    cdef int M = X.shape[0] # M rows
    cdef int N = X.shape[1] # N columns

    # simulate doublet indices
    np.random.seed(random_state)
    cdef Py_ssize_t n_sim = <int>(M * sim_doublet_ratio)
    cdef int[:, :] doublet_indices = np.random.randint(0, M, size=(n_sim, 2), dtype = np.int32)

    # if sparse
    cdef const int[:] data
    cdef const int[:] indices
    cdef const long[:] indptr

    cdef int[:] out_data
    cdef int[:] out_indices
    cdef long[:] out_indptr

    # if not sparse
    cdef const int[:, :] data_array
    cdef int[:, :] out_array

    # Generate new count matrix
    if issparse(X):
        Xdata = X.data
        if Xdata.dtype != np.int32:
            Xdata = Xdata.astype(np.int32)

        data = Xdata
        indices = X.indices
        indptr = X.indptr.astype(np.int64)

        size = 0
        for i in range(n_sim):
            u = doublet_indices[i, 0]
            v = doublet_indices[i, 1]
            size += (indptr[u + 1] - indptr[u]) + (indptr[v + 1] - indptr[v])

        out_data_buffer = np.zeros(size, dtype = np.int32)
        out_indices_buffer = np.zeros(size, dtype = np.int32)
        out_indptr_buffer = np.zeros(n_sim + 1, dtype = np.int64)
        
        out_data = out_data_buffer
        out_indices = out_indices_buffer
        out_indptr = out_indptr_buffer

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

        results = csr_matrix((out_data_buffer, out_indices_buffer, out_indptr_buffer), shape = (n_sim, N), copy = False)
    else:
        results = np.zeros((n_sim, N), dtype = np.int32)
        data_array = X
        out_array = results

        for i in range(n_sim):
            u = doublet_indices[i, 0]
            v = doublet_indices[i, 1]
            for j in range(N):
                out_array[i, j] = data_array[u, j] + data_array[v, j]

    return results, np.asarray(doublet_indices)
