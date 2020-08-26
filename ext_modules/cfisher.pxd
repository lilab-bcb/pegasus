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

cdef double pvalue(int a_true, int a_false, int b_true, int b_false) nogil

cdef void cfisher(int s, const int* a_true, const int* a_false, const int* b_true, const int* b_false, float* oddsratio, float* pval) nogil
