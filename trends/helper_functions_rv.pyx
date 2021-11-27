# cython: language_level=3, boundscheck=False, cdivision=True, wraparound=False
# cython: binding=True

import numpy as np
cimport numpy as np
from c_kepler._kepler import kepler_single

from astropy.time import Time
from tqdm import tqdm
import cython
from libc.math cimport sin, cos, tan, atan, sqrt, log


cdef float two_pi, math_e, G, auday2ms, hip_beginning

two_pi = 6.283185307179586
math_e = 2.718281828459045

# G in AU, M_Jup, day units.
G = 2.824760877012879e-07 # (c.G.cgs*(1/c.au.cgs)**3 * (c.M_jup.cgs) * (24*3600)**2).value

# Converts AU/day to m/s
auday2ms = 1731456.8368055555 # c.au.si.value/(24*3600)

# Just need the "zero time" to evolve mean anomaly into the rv_epoch
hip_beginning = Time(1989.85, format='decimalyear').jd

def rv_list(double [:] a_list, double [:] m_list, double [:] e_list, 
            double [:] i_list, double [:] om_list, double [:] M_anom_0_list,
            double [:] per_list, double m_star, double rv_epoch, 
            double gdot, double gdot_err, double gddot, double gddot_err):
    
    cdef int num_points, j
    cdef double a, m, e, i, om, M_anom_0, per, log_lik
    num_points = a_list.shape[0]
               
    cdef np.ndarray[double, ndim=1] lik_list = np.ndarray(shape=(num_points,),
                                                            dtype=np.float64)
    
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
    
    Returns:
        log_likelihood_total (float): Likelihood of gdot AND gddot given the model.
    """
    cdef double E, gdot_model, gddot_model
    cdef double log_likelihood_gdot, log_likelihood_gddot, log_likelihood_total
    
    #per = hlp.P(a, m, m_star)
    E = M_2_evolvedE(M_anom_0, per, e, rv_epoch)

    gdot_model, gddot_model = gamma(a, m, e, i, om, E, per, m_star)

    # Log of the prefactor minus the log of the exponential
    log_likelihood_gdot  = log(1/(sqrt(two_pi)*gdot_err))\
                         - (gdot-gdot_model)**2/(2*gdot_err**2)
                         
    log_likelihood_gddot = log(1/(sqrt(two_pi)*gddot_err))\
                         - (gddot-gddot_model)**2/(2*gddot_err**2)
                 
    log_likelihood_total = log_likelihood_gdot + log_likelihood_gddot

    return log_likelihood_total


# It seems gamma() needs to be a cdef function, otherwise it returns nans
cpdef (double, double) gamma(double a, double m, double e, 
                             double i, double om, double E, 
                             double per, double m_star):
    """
    Given an orbital model, calculate gdot and gddot.
    
    Arguments:
        a (float, AU): Semi-major axis
        m (float, M_jup): Companion mass
        e (float): Orbital eccentricity
        i (float, radians): Orbital inclination
        om (float, radians): Argument of periastron
        E (float, radians): Eccentric anomaly
        per (float, days): Orbital period
        m_star (float, M_jup): Host star mass
    
    Returns:
        gdot (float, m/s/day): Linear trend term
        gddot (float, m/s/day/day): Quadratic curvature term
    """

    cdef double     e_term, sqrt_eterm,\
                    sqrt_e_sq_term, cos_E, sin_E,\
                    tan_Eovr2, nu, nu_dot, nu_ddot,\
                    cos_nu_om, sin_nu_om, sin_i,\
                    pre_fac, gamma_dot, gamma_ddot

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
    pre_fac = sqrt(G)/sqrt_e_sq_term * m*sin_i/sqrt((m+m_star)*(a)) * auday2ms # AU/day ==> m/s

    gamma_dot = -pre_fac*nu_dot*sin_nu_om # m/s/day
    gamma_ddot = -pre_fac*(nu_dot**2*cos_nu_om + nu_ddot*sin_nu_om) # m/s/day/day

    return gamma_dot, gamma_ddot

    
cpdef M_2_evolvedE(double M0, double per, double e, double rv_epoch):
    """
    Takes a mean anomaly at the beginning of the Hip mission, 
    evolves it to the RV epoch, and converts it to eccentric anomaly.

    M0 is the mean anomaly in radians.
    per is the period in days.
    e is the eccentricity.
    rv_epoch is the bjd corresponding to the ~midpoint of the RV baseline, 
    where gdot and gddot are measured.
    """
    # Fewer python references with no declarations. Not sure why.
    
    M_evolved = ((two_pi/per)*(rv_epoch - hip_beginning) + M0)%two_pi

    E = kepler_single(M_evolved, e)

    return E
    