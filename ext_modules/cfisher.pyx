"""
Codes copied and slighly modified from https://github.com/brentp/fishers_exact_test

Original license is BSD 3-Clause (https://github.com/brentp/fishers_exact_test/blob/master/LICENSE)
"""

"""
Cython Fisher's exact test:
Fisher's exact test is a statistical significance test used in the
analysis of contingency tables where sample sizes are small.
Function lngamma(), lncombination(), hypergeometric_probability(),
were originally written by Oyvind Langsrud:
Oyvind Langsrud
Copyright (C) : All right reserved.
Contact Oyvind Langsrud for permission.
Adapted to Cython version by:
Haibao Tang, Brent Pedersen
"""

# cython: language_level=3

import numpy as np

cimport cython

cdef extern from "math.h":
    double log(double) nogil
    double exp(double) nogil
    double lgamma(double) nogil


cdef inline double _naive_lnfactorial(int n) nogil:
    cdef double acc = 0.0
    cdef int i
    for i in range(2, n + 1):
        acc += log(i)
    return acc

# Tabulated ln n! for n \in [0, 1023]
cdef double[1024] _lnfactorials1
cdef int i
for i in range(1024):
    _lnfactorials1[i] = _naive_lnfactorial(i)


# Logarithm of n! with algorithmic approximation
@cython.boundscheck(False)
cdef inline double lnfactorial(int n) nogil:
    return _lnfactorials1[n] if n < 1024 else lgamma(n + 1)


# Logarithm of the number of combinations of 'n' objects taken 'p' at a time
cdef inline double lncombination(int n, int p) nogil:
    return lnfactorial(n) - lnfactorial(p) - lnfactorial(n - p)


# Compute the hypergeometric probability, or probability that a list of
# 'n' objects should contain 'x' ones with a particular property when the
# list has been selected randomly without replacement from a set of 'N'
# objects in which 'K' exhibit the same property
cdef inline double hypergeometric_probability(int x, int n, int K, int N) nogil:
    return exp(lncombination(K, x)
               + lncombination(N - K, n - x)
               - lncombination(N, n))


cdef double pvalue(int a_true, int a_false, int b_true, int b_false) nogil:
    # Convert the a/b groups to study vs population.
    cdef int k = a_true
    cdef int n = a_false + a_true  # total in study.
    cdef int K = a_true + b_true
    cdef int N = K + a_false + b_false

    cdef int lm = max(0, n - (N - K))
    cdef int um = min(n, K)
    if lm == um:
        return 1.0

    cdef double epsilon = 1e-6
    cdef double cutoff = hypergeometric_probability(k, n, K, N)
    cdef double two_tail = 0
    cdef int i
    cdef double p

    for x in range(lm, um + 1):
        p = hypergeometric_probability(x, n, K, N)
        if p <= cutoff + epsilon:
            two_tail += p

    return min(two_tail, 1.0)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef void cfisher(int s, const int* a_true, const int* a_false, const int* b_true, const int* b_false, float* oddsratio, float* pval) nogil:
    cdef int i
    for i in range(s):
        # oddsratio
        if b_true[i] == 0 or a_false[i] == 0:
            oddsratio[i] = 0.0 if a_true[i] == 0 or b_false[i] == 0 else 1e30
        else:
            oddsratio[i] = (a_true[i] * 1.0 / a_false[i]) * (b_false[i] * 1.0 / b_true[i])
        # two-sided p value
        pval[i] = pvalue(a_true[i], a_false[i], b_true[i], b_false[i])


def fisher_exact(const int[:] a_true, const int[:] a_false, const int[:] b_true, const int[:] b_false):
    cdef int s = a_true.size

    oddsratio = np.zeros(s, dtype = np.float32)
    pval = np.zeros(s, dtype = np.float32)

    cdef float[:] orview = oddsratio
    cdef float[:] pview = pval

    cfisher(s, &a_true[0], &a_false[0], &b_true[0], &b_false[0], &orview[0], &pview[0])

    return oddsratio, pval
