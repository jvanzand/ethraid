import numpy as np
cimport numpy as np
from libc.math cimport sin, cos, tan


#import scipy as sp
#import astropy.constants as c
#import astropy.units as u
import cython

#import radvel as rv

#from c_kepler import _kepler as ck

#from kern_profiler_dummy import *

## Constants ##

cdef float pi, G, M_sun, M_jup, au

pi = 3.141592653589793
G =  6.674299999999999e-08
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
au = 14959787070000.0


cpdef gamma(double a, double Mp, double per, double e, double i, double om, double E):
    
    cdef double Mp_units, a_units, nu, E_dot, nu_dot, prefac, gd_t1, gd_t2, gamma_dot, gd_t1_dot, gd_t2_dot, gdd_t1, gdd_t2, gamma_ddot
    
    Mp_units = Mp*M_jup
    a_units = a*au
    
    nu = 2*np.arctan(((1+e)/(1-e))**0.5*tan(E/2))

    # Differentiate Kepler's equation in time to get E_dot
    # Note that E_dot has units of (1/per), where [per] is days. Therefore [gamma_ddot] = m/s/d^2
    E_dot = (2*pi/per)/(1-e*cos(E))
    nu_dot = (1+tan(nu/2)**2)**-1 * ((1+e)/(1-e))**0.5 * cos(E/2)**-2 * E_dot

    # Convert prefac units from cm/s^2 to m/s/day
    # Negative just depends on choice of reference direction. I am being consistent with radvel rv_drive function.
    prefac = -(Mp_units*G*sin(i))/(a_units**2*(1-e)) * (1/100) * (24*3600)


    gd_t1 = (1+cos(nu))/(1+cos(E))
    gd_t2 = sin(nu+om)/(1-e*cos(E))


    gamma_dot = prefac*gd_t1*gd_t2

    gd_t1_dot = ((1+cos(nu))*sin(E) * E_dot - (1+cos(E))*sin(nu)*nu_dot) / (1+cos(E))**2
    gd_t2_dot = ((1-e*cos(E))*cos(nu+om) * nu_dot - sin(nu+om)*e*sin(E)*E_dot) / (1-e*cos(E))**2


    gdd_t1 = gd_t2 * gd_t1_dot
    gdd_t2 = gd_t1 * gd_t2_dot

    gamma_ddot = prefac*(gdd_t1+gdd_t2)

    return gamma_dot, gamma_ddot

cpdef T_progression(double T_anom_0, double e, double [:] E_prog, int t_num):

    cdef int i
    
    cdef np.ndarray[double, ndim=1] T_prog = \
        np.ndarray(shape=(t_num,), dtype=np.float64)
    
    for i in range(t_num):
        T_prog[i] = T_anom_0 + 2*np.arctan( ((1+e)/(1-e))**0.5 * tan(E_prog[i]/2))


    return T_prog








