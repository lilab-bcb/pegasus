# cython: language_level=3

import numpy as np
import scipy.stats as ss
cimport cython

from libc.math cimport sqrt, log2, M_LOG2E, fabs
# from libc.stdio cimport printf

ctypedef fused indices_type:
    int
    long

ctypedef fused indptr_type:
    int
    long

# Note that for now csr_to_csc or csr_to_csc_cond return long for indices and indptr; Once we update Cython to 3.0, we will use const fused memoryview instead!!!

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple csr_to_csc(const float[:] input_data, indices_type[:] input_indices, indptr_type[:] input_indptr, int M, int N, const long[:] ords):
    """ This routine in addition group cells by clusters (ords)"""
    cdef Py_ssize_t i, j, pos, col

    output_indptr = np.zeros(N+1, dtype = np.int64)
    cdef long[:] indptr = output_indptr
    cdef long[:] counter = np.zeros(N, dtype = np.int64)

    for i in range(input_indices.size):
        if input_data[i] != 0.0: # in case there are extra 0s in the sparse matrix
            indptr[input_indices[i]+1] += 1
    for i in range(N):
        counter[i] = indptr[i]
        indptr[i + 1] += indptr[i]

    output_data = np.zeros(indptr[N], dtype = np.float32)
    output_indices = np.zeros(indptr[N], dtype = np.int64)
    cdef float[:] data = output_data
    cdef long[:] indices = output_indices

    for i in range(M):
        for j in range(input_indptr[ords[i]], input_indptr[ords[i]+1]):
            if input_data[j] != 0.0:
                col = input_indices[j]
                pos = counter[col]
                data[pos] = input_data[j]
                indices[pos] = i
                counter[col] += 1

    return output_data, output_indices, output_indptr


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef tuple csr_to_csc_cond(const float[:] input_data, indices_type[:] input_indices, indptr_type[:] input_indptr, int N, const long[:] ords, const long[:] cumsum):
    cdef Py_ssize_t i, j, k, pos, col
    cdef Py_ssize_t n, start, end, fr, to

    n = cumsum.size

    output_data_list = []
    output_indices_list = []
    output_indptrs = np.zeros((n, N+1), dtype = np.int64)

    cdef float[:] data
    cdef long[:] indices
    cdef long[:] indptr
    cdef long[:] counter = np.zeros(N, dtype = np.int64)

    start = 0
    for i in range(n):
        end = cumsum[i]

        indptr = output_indptrs[i]

        for j in range(start, end):
            fr = input_indptr[ords[j]]
            to = input_indptr[ords[j] + 1]
            for k in range(fr, to):
                if input_data[k] != 0.0:
                    indptr[input_indices[k] + 1] += 1

        for j in range(N):
            counter[j] = indptr[j]
            indptr[j + 1] += indptr[j]

        output_data_list.append(np.zeros(indptr[N], dtype = np.float32))
        output_indices_list.append(np.zeros(indptr[N], dtype = np.int64))

        data = output_data_list[i]
        indices = output_indices_list[i]
        for j in range(start, end):
            for k in range(input_indptr[ords[j]], input_indptr[ords[j] + 1]):
                if input_data[k] != 0.0:
                    col = input_indices[k]
                    pos = counter[col]
                    data[pos] = input_data[k]
                    indices[pos] = j - start
                    counter[col] += 1

        start = end

    return output_data_list, output_indices_list, output_indptrs


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void partition_indices(const long[:] indices, const long[:] cumsum, long[:] posarr):
    cdef Py_ssize_t i, j, s

    posarr[0] = 0
    i = 0
    s = indices.size

    for j in range(cumsum.size):
        while i < s and indices[i] < cumsum[j]:
            i += 1
        posarr[j + 1] = i


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef tuple calc_mwu(int start_pos, int end_pos, const float[:] data, const long[:] indices, const long[:] indptr, const long[:] n1arr, const long[:] n2arr, const long[:] cumsum, int first_j, int second_j, bint verbose):
    """ Run Mann-Whitney U test for all clusters in cluster_labels, focusing only on genes in [start_pos, end_pos)
        P values are calculated based on normal distribution approximation, accurate when each cluster has > 20 samples
        This function assumes the cells are grouped by cluster ids
    """
    cdef Py_ssize_t ngene, ncluster, nsample
    cdef Py_ssize_t i, j, pos_i
    cdef Py_ssize_t n_nonzero, n_zero
    cdef double sd_factor, avg_zero_rank
    cdef double R1, U1, U2, mean_U, sd_U
    
    ngene = end_pos - start_pos
    ncluster = cumsum.size
    nsample = cumsum[ncluster - 1]

    U_stats_np = np.zeros((ngene, ncluster), dtype = np.float64)
    pvals_np = np.ones((ngene, ncluster), dtype = np.float32)
    aurocs_np = np.full((ngene, ncluster), 0.5, dtype = np.float32)
    # memoryviews
    cdef double[:, :] U_stats = U_stats_np
    cdef float[:, :] pvals = pvals_np
    cdef float[:, :] aurocs = aurocs_np

    # one-vs-rest at cluster level
    n1Rs_np = np.zeros(ncluster, dtype = np.float64)
    n1n2_np = np.zeros(ncluster, dtype = np.float64)
    U_default_np = np.zeros(ncluster, dtype = np.float64)
    # memoryviews
    cdef double[:] n1Rs = n1Rs_np
    cdef double[:] n1n2 = n1n2_np
    cdef double[:] U_default = U_default_np

    cdef long[:] posarr = np.zeros(ncluster + 1, dtype = np.int64)

    cdef Py_ssize_t buffer_size = ngene * ncluster
    cdef Py_ssize_t buffer_pos = 0
    cdef double[:] zscores = np.zeros(buffer_size, dtype = np.float64)

    for j in range(ncluster):
        n1Rs[j] = (<double>(n1arr[j] + 1)) * n1arr[j] / 2.0
        n1n2[j] = (<double>n1arr[j]) * n2arr[j]
        U_default[j] = n1n2[j] / 2.0

    for i in range(start_pos, end_pos):
        pos_i = i - start_pos
        n_nonzero = indptr[i + 1] - indptr[i]
        n_zero = nsample - n_nonzero
        if n_nonzero == 0:
            U_stats_np[pos_i] = U_default_np
        else:
            ranks = ss.rankdata(data[indptr[i]:indptr[i + 1]]) + n_zero # np.float64

            # (n+1) - \sum_{i=1}^{k} \frac{t_i^3 - t_i}{n(n-1)}
            _, ties = np.unique(ranks, return_counts = True)
            if n_zero > 0:
                ties = np.concatenate(([n_zero], ties))
            sd_factor = 0.0 # sd factor for tie correction
            if nsample > 1:
                sd_factor = sqrt(((nsample + 1) - (ties**3.0 - ties).sum() / (nsample - 1.0) / nsample) / 12.0)

            avg_zero_rank = (n_zero + 1.0) / 2.0
            partition_indices(indices[indptr[i]: indptr[i + 1]], cumsum, posarr)
            for j in range(ncluster):
                if n1arr[j] == 0 or n2arr[j] == 0:
                    U_stats[pos_i, j] = U_default[j]
                elif j == second_j:
                    U_stats[pos_i, j] = n1n2[first_j] - U_stats[pos_i, first_j]
                    aurocs[pos_i, j] = 1.0 - aurocs[pos_i, first_j]
                else:
                    R1 = ranks[posarr[j]:posarr[j+1]].sum() + (n1arr[j] - (posarr[j+1] - posarr[j])) * avg_zero_rank
                    U1 = R1 - n1Rs[j]
                    U2 = n1n2[j] - U1
                    mean_U = U_default[j] + 0.5 # continuity correction
                    sd_U = sqrt(n1n2[j]) * sd_factor
                    z = (max(U1, U2) - mean_U) / sd_U

                    zscores[buffer_pos] = z
                    buffer_pos += 1

                    U_stats[pos_i, j] = U1
                    aurocs[pos_i, j] = U1 / n1n2[j]

    cdef double[:] buffer_pvals = ss.norm.sf(zscores[0:buffer_pos]) * 2.0
    buffer_pos = 0
    for i in range(start_pos, end_pos):
        if indptr[i + 1] > indptr[i]:
            for j in range(ncluster):
                if n1arr[j] > 0 and n2arr[j] > 0 and j != second_j:
                    pvals[i - start_pos, j] = buffer_pvals[buffer_pos]
                    buffer_pos += 1

    if second_j > 0:
        pvals[:, second_j] = pvals[:, first_j]

    if verbose:
        print(f"calc_mwu finished for genes in [{start_pos}, {end_pos}).")

    return U_stats_np, pvals_np, aurocs_np



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cpdef list calc_stat(const float[:] data, const long[:] indices, const long[:] indptr, const long[:] n1arr, const long[:] n2arr, const long[:] cumsum, int first_j, int second_j, bint t_test, bint fisher_test, bint verbose):
#     """ Collect sufficient statistcs and optionally run Welch's T test and Fisher's exact test
#     """
    cdef Py_ssize_t ngene, ncluster, buffer_size
    cdef Py_ssize_t i, j, k

    ngene = indptr.size - 1
    ncluster = cumsum.size
    buffer_size = ngene * ncluster

    sstats_np = np.zeros((6, ngene, ncluster), dtype = np.float32) # log2Mean, log2Mean_other, log2FC, percentage, percentage_other, percentage_fold_change
    cdef float[:,:,:] sstats = sstats_np


    #sufficient statistics
    cdef double tot_sum
    cdef Py_ssize_t tot_nzero
    cdef double[:] sums = np.zeros(ncluster, dtype = np.float64)
    cdef double[:] mean1 = np.zeros(ncluster, dtype = np.float64)
    cdef double[:] mean2 = np.zeros(ncluster, dtype = np.float64)
    cdef long[:] nonzero = np.zeros(ncluster, dtype = np.int64)
    cdef long[:] posarr = np.zeros(ncluster + 1, dtype = np.int64)
    #Welch's t-test
    cdef double M_LOG2E2 = M_LOG2E * M_LOG2E
    cdef double tot_sum2, s1sqr, s2sqr, var_est
    cdef double[:] sum2s = np.zeros(ncluster, dtype = np.float64)
    
    cdef float[:,:,:] t_arr # t_tstat, t_pval
    cdef Py_ssize_t ttest_pos = 0
    cdef double[:] tscore, upsilon
    cdef double[:] buffer_pvals
    #Fisher's exact test
    cdef int[:,:] a_true, a_false, b_true, b_false


    if t_test:
        t_arr_np = np.zeros((2, ngene, ncluster), dtype = np.float32)
        t_arr = t_arr_np
        tscore = np.zeros(buffer_size, dtype = np.float64)
        upsilon = np.zeros(buffer_size, dtype = np.float64)

    if fisher_test:
        a_true_np = np.zeros((ncluster, ngene), dtype = np.int32)
        a_false_np = np.zeros((ncluster, ngene), dtype = np.int32)
        b_true_np = np.zeros((ncluster, ngene), dtype = np.int32)
        b_false_np = np.zeros((ncluster, ngene), dtype = np.int32)
        a_true = a_true_np
        a_false = a_false_np
        b_true = b_true_np
        b_false = b_false_np


    for i in range(ngene):
        tot_sum = 0.0
        tot_nzero = 0
        partition_indices(indices[indptr[i]:indptr[i+1]], cumsum, posarr)

        for j in range(ncluster):
            sums[j] = 0.0
            nonzero[j] = 0
            for k in range(indptr[i] + posarr[j], indptr[i] + posarr[j+1]):
                sums[j] += data[k]
                nonzero[j] += 1
            sums[j] *= M_LOG2E
            tot_sum += sums[j]
            tot_nzero += nonzero[j]

        for j in range(ncluster):
            if j == second_j:
                sstats[0, i, j] = sstats[1, i, first_j]
                sstats[1, i, j] = sstats[0, i, first_j]
                sstats[2, i, j] = -sstats[2, i, first_j]
                sstats[3, i, j] = sstats[4, i, first_j]
                sstats[4, i, j] = sstats[3, i, first_j]
                sstats[5, i, j] = 0 if sstats[5, i, first_j] == 1e30 else (1.0 / sstats[5, i, first_j] if sstats[5, i, first_j] > 0.0 else 1e30)
            else:
                if n1arr[j] > 0:
                    mean1[j] = sums[j] / n1arr[j]
                    sstats[0, i, j] = mean1[j]
                    sstats[3, i, j] = nonzero[j] * 100.0 / n1arr[j]
                if n2arr[j] > 0:
                    mean2[j] = max(tot_sum - sums[j], 0.0) / n2arr[j]
                    sstats[1, i, j] = mean2[j]
                    sstats[4, i, j] = (tot_nzero - nonzero[j]) * 100.0 / n2arr[j]
                sstats[2, i, j] = sstats[0, i, j] - sstats[1, i, j]
                sstats[5, i, j] = sstats[3, i, j] / sstats[4, i, j] if sstats[4, i, j] > 0.0 else 1e30

        if t_test:
            tot_sum2 = 0.0
            for j in range(ncluster):
                sum2s[j] = 0.0
                for k in range(indptr[i] + posarr[j], indptr[i] + posarr[j+1]):
                    sum2s[j] += <double>(data[k]) ** 2.0
                sum2s[j] *= M_LOG2E2
                tot_sum2 += sum2s[j]

            for j in range(ncluster):
                s1sqr = s2sqr = 0.0
                t_arr[1, i, j] = 1.0
                if n1arr[j] > 1 and n2arr[j] > 1 and j != second_j:
                    s1sqr = max(sum2s[j] - n1arr[j] * (mean1[j] ** 2.0), 0.0) / (n1arr[j] - 1.0)
                    s2sqr = max((tot_sum2 - sum2s[j]) - n2arr[j] * (mean2[j] ** 2.0), 0.0) / (n2arr[j] - 1.0)

                    var_est = s1sqr / n1arr[j] + s2sqr / n2arr[j]
                    if var_est > 0.0:
                        tscore[ttest_pos] = (mean1[j] - mean2[j]) / sqrt(var_est)
                        t_arr[0, i, j] = tscore[ttest_pos]
                        tscore[ttest_pos] = fabs(tscore[ttest_pos])
                        upsilon[ttest_pos] = (var_est ** 2.0) / ((s1sqr / n1arr[j]) ** 2.0 / (n1arr[j] - 1.0) + (s2sqr / n2arr[j]) ** 2.0 / (n2arr[j] - 1.0))
                    else:
                        upsilon[ttest_pos] = 1.0
                    ttest_pos += 1

        if fisher_test:
            for j in range(ncluster):
                if j != second_j:
                    a_true[j, i] = nonzero[j]
                    a_false[j, i] = n1arr[j] - a_true[j, i]
                    b_true[j, i] = tot_nzero - nonzero[j]
                    b_false[j, i] = n2arr[j] - b_true[j, i]

    results = [sstats_np]

    if t_test:
        buffer_pvals = ss.t.sf(tscore[0:ttest_pos], upsilon[0:ttest_pos]) * 2.0 # two-sided
        ttest_pos = 0
        for i in range(ngene):
            for j in range(ncluster):
                if j == second_j:
                    t_arr[0, i, j] = -t_arr[0, i, first_j]
                    t_arr[1, i, j] = t_arr[1, i, first_j]
                elif n1arr[j] > 1 and n2arr[j] > 1 and j != second_j:
                    t_arr[1, i, j] = buffer_pvals[ttest_pos]
                    ttest_pos += 1
        results.append(t_arr_np)

    if fisher_test:
        results.append([a_true_np, a_false_np, b_true_np, b_false_np])

    return results
