#import pyximport
#pyximport.install()

import numpy as np

cdef class best_class:
    cdef float height, width
	
    def __init__(self, h, w):
        self.height=h
        self.width=w
		
    def attributes(self):
        print('I am a tree and my height is ', self.height, self.width)
        
    def dense_loop(self, long long n):
        cdef unsigned int i
        my_list = []
        for i in range(n):
            if i%17 == 0:
                my_list.append(i)
        return my_list

    def dense_loopp(self, n):
        my_list = []
        for i in range(n):
            if i%17 == 0:
                my_list.append(i)
        return my_list
		
		
class python_class(best_class):
    pass
