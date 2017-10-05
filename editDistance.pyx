# Cython module for edit_distance function
# Define types with c/cp-def
# Convert Python parameters to those of C (test -> a)
# Read through length (l) of string and count mismatches (m)

cpdef int edit_distance(test, ref):
    cdef char * a = test
    cdef char * b = ref
    cdef int k, l, m

    m = 0
    l = len(a)

    for k from 0 <= k < l:
        if a[k] != b[k]:
            m += 1

    return m