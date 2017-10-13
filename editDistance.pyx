# Cython module for edit_distance function
# Define types with c/cp-def and convert Python parameters to those of C (test -> a)
# To increase performance, convert Python string to C bytes. Some type conversion required for Python 3
# Read through length (l) of string and count mismatches (m)

cpdef int edit_distance(test, ref):
    cdef bytes a = <bytes>test
    cdef bytes b = <bytes>ref
    cdef int k, l, m

    m = 0
    l = len(a)

    for k from 0 <= k < l:
        if a[k] != b[k]:
            m += 1

    return m