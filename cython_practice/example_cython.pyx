from libc.stdio cimport printf

cpdef void test(int x):
    cdef int y = 0
    cdef int i
    
    for i in range(x):
        y += i
    #print(y)
        printf('%i ', y)
    return