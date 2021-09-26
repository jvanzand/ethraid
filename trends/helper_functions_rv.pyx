# cython: language_level=3, boundscheck=False, cdivision=True, wraparound=False
# cython: binding=True

#import os
#import sys
#path = os.getcwd()
#sys.path.append(path+'/trends') # To import c_kepler
import numpy as np
cimport numpy as np
from c_kepler._kepler import kepler_single

from astropy.time import Time
from tqdm import tqdm
import cython
from libc.math cimport sin, cos, tan, atan, sqrt, log

import helper_functions_rv as hlp_rv

##########################################
#### Kepler solver for one M and one e
# Wrapping kepler(M,e) a simple function that takes two doubles as
# arguments and returns a double
#cdef extern from "../c_kepler/kepler.c":
#    double kepler(double M, double e)
#    double rv_drive(double t, double per, double tp, double e, double cosom, double sinom, double k )
##########################################

cdef float pi, two_pi, math_e, G, M_sun, M_jup, au, pc_in_cm, hip_beginning

pi = 3.141592653589793
two_pi = 6.283185307179586
math_e = 2.718281828459045
G = 6.674299999999999e-08
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
au = 14959787070000.0
pc_in_cm = 3.086e18

# Just need the "zero time" to evolve mean anomaly into the rv_epoch
hip_beginning = Time(1989.85, format='decimalyear').jd

def rv_list(double [:] a_list, double [:] m_list, double [:] e_list, 
               double [:] i_list, double [:] om_list, double [:] M_anom_0_list,
               double [:] per_list, double m_star, double rv_epoch, 
               double gdot, double gdot_err, double gddot, double gddot_err):
    
    cdef int num_points, j
    cdef double a, m, e, i, om, M_anom_0, per, log_lik
    num_points = a_list.shape[0]
               
    cdef np.ndarray[double, ndim=1] lik_list = np.ndarray(shape=(num_points,), dtype=np.float64)
    
    print('Running RV models')
    for j in tqdm(range(num_points)):
       
        a  =  a_list[j]
        m  =  m_list[j]
        e  =  e_list[j]
        i  =  i_list[j]
        om = om_list[j]
        M_anom_0  =  M_anom_0_list[j]
        per = per_list[j]
        
        log_lik = log_lik_gamma(a, m, e, i, om, M_anom_0,
                                per, m_star, rv_epoch,
                                gdot, gdot_err, gddot, gddot_err)
                              
        lik_list[j] = math_e**log_lik
    return lik_list

def log_lik_gamma(double a, double m, double e, double i, double om, double M_anom_0, 
                  double per, double m_star, double rv_epoch,
                  double gdot, double gdot_err, double gddot, double gddot_err):
    """
    Compute the log-likelihood of a state (set of a, Mp, e, i, om, and M_anom_0)
    given the RV data (true gammas and their uncertainties).
    """
    
    #cdef double a, m, e, i, om, M_anom_0 
    #cdef double gdot_data, gdot_err, gddot_data, gddot_err
    cdef double E, gdot_model, gddot_model
    cdef double log_likelihood_gdot, log_likelihood_gddot, log_likelihood_total
    
    #per = hlp.P(a, m, m_star)
    E = hlp_rv.M_2_evolvedE(M_anom_0, per, e, rv_epoch)

    gdot_model, gddot_model = gamma(a, m, e, i, om, E, per, m_star)

    # Log of the prefactor minus the log of the exponential
    log_likelihood_gdot  = log(1/(sqrt(two_pi)*gdot_err))\
                         - (gdot-gdot_model)**2/(2*gdot_err**2)
                         
    log_likelihood_gddot = log(1/(sqrt(two_pi)*gddot_err))\
                         - (gddot-gddot_model)**2/(2*gddot_err**2)
                 
    log_likelihood_total = log_likelihood_gdot + log_likelihood_gddot

    return log_likelihood_total


# It seems gamma() needs to be a cdef function, otherwise it returns nans
# Testing the above comment by making it a cpdef so I can use it in log_likelihood.py
#@profile
cpdef (double, double) gamma(double a, double m, double e, 
                             double i, double om, double E, 
                             double per, double m_star):

    cdef double     m_g, m_star_g, a_cm, e_term, sqrt_eterm,\
                    sqrt_e_sq_term, cos_E, sin_E,\
                    tan_Eovr2, nu, nu_dot, nu_ddot,\
                    cos_nu_om, sin_nu_om, sin_i,\
                    pre_fac, gamma_dot, gamma_ddot

    #per_sec = per*86400 # 24*3600 to convert days ==> seconds
    m_g = m*M_jup
    m_star_g = m_star*M_sun

    a_cm = a*au

    e_term = (1+e)/(1-e)
    sqrt_eterm = sqrt(e_term)
    sqrt_e_sq_term = sqrt(1-e*e)

    cos_E = cos(E)
    sin_E = sin(E)
    tan_Eovr2 = sin_E/(1+cos_E)

    nu = 2*atan(sqrt_eterm*tan_Eovr2)

    # nu derivatives use days (not seconds) to give gdot/gddot correct units 
    nu_dot = two_pi*sqrt_e_sq_term/(per*(1-e*cos_E)**2) # Units of day^-1
    nu_ddot = -nu_dot**2 * 2*e*sin_E/sqrt_e_sq_term # # Units of day^-2

    cos_nu_om = cos(nu+om)
    sin_nu_om = sin(nu+om)
    sin_i = sin(i)

    # Fischer (analytic)
    pre_fac = sqrt(G)/sqrt_e_sq_term * m_g*sin_i/sqrt((m_g+m_star_g)*(a_cm)) * 1/100 # cm/s ==> m/s

    gamma_dot = -pre_fac*nu_dot*sin_nu_om # m/s/d
    gamma_ddot = -pre_fac*(nu_dot**2*cos_nu_om + nu_ddot*sin_nu_om) # m/s/d/d

    return gamma_dot, gamma_ddot

    
cpdef M_2_evolvedE(double M0, double per, double e, double rv_epoch):
    """
    Takes a mean anomaly at the beginning of the Hip mission, 
    evolves it to the RV epoch, and converts it to eccentric anomaly.

    M0 is the mean anomaly in radians.
    per is the period in days
    e is the eccentricity
    rv_epoch is the bjd corresponding to the ~midpoint of the RV baseline. 
    It is where gdot and gddot are measured
    """
    # Fewer python references with no declarations
    
    M_evolved = ((two_pi/per)*(rv_epoch - hip_beginning) + M0)%two_pi

    E = kepler_single(M_evolved, e)

    return E
    