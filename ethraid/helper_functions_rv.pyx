import numpy as np
from tqdm import tqdm
import cython
cimport numpy as np
from libc.math cimport sin, cos, tan, atan, sqrt, log

from ethraid import hip_times
from ethraid.compiled._kepler import kepler_single


cdef double two_pi, math_e, G, auday2ms, hip_beginning

two_pi = 6.283185307179586
math_e = 2.718281828459045

# G in {AU, M_Jup, day} units.
G = 2.824760877012879e-07 # (c.G.cgs*(1/c.au.cgs)**3 * (c.M_jup.cgs) * (24*3600)**2).value

# Converts AU/day to m/s
auday2ms = 1731456.8368055555 # c.au.si.value/(24*3600)

# Just need the "zero time" to evolve mean anomaly into the rv_epoch
hip_beginning = hip_times[0]

def rv_list(double [:] a_list, double [:] m_list, double [:] e_list, 
            double [:] i_list, double [:] om_list, double [:] M_anom_0_list,
            double [:] per_list, double m_star, double rv_epoch, 
            double gdot, double gdot_err, double gddot, double gddot_err):
    
    """
    Calculates the likelihood of the RV data (measured gdot, gddot,
    and their uncertainties) conditioned on each of a list of orbital 
    models, resulting in a list of likelihoods. All input lists must 
    have the same length, and the output lik_list has that length as well.
    
    Arguments:
        a_list (list of floats, AU): Semi-major axis
        m_list (list of floats, M_jup): Companion mass
        e_list (list of floats): Orbital eccentricity
        i_list (list of floats, radians): Orbital inclination
        om_list (list of floats, radians): Argument of periastron
        M_anom_0_list (list of floats, radians): 
        per_list (list of floats, days): Orbital period
        m_star (float, M_jup): Host star mass
        rv_epoch (float, days): BJD at which the model RV trend
                                and curv are calculated.
        gdot (float, m/s/day): Measured linear trend term
        gdot (float, m/s/day): Error on gdot
        gddot (float, m/s/day/day): Measured quadratic curvature 
                                    term
        gddot_err (float, m/s/day/day): Error on gddot
    """
    
    cdef int num_points, j
    cdef double a, m, e, i, om, M_anom_0, per, log_lik
    num_points = a_list.shape[0]
               
    cdef np.ndarray[double, ndim=1] lik_list = np.ndarray(shape=(num_points,),
                                                            dtype=np.float64)
    if gdot_err==0 or gddot_err==0:
        raise Exception('Errors cannot be 0')
        
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
    Computes the log-likelihood of the RV data conditioned on
    a model (set of a, Mp, e, i, om, and M_anom_0).
    
    Arguments:
        a (float, AU): Semi-major axis
        m (float, M_jup): Companion mass
        e (float): Orbital eccentricity
        i (float, radians): Orbital inclination
        om (float, radians): Argument of periastron
        M_anom_0 (float, radians): 
        per (float, days): Orbital period
        m_star (float, M_jup): Host star mass
        rv_epoch (float, days): BJD at which the model RV trend
                                and curv are calculated.
        gdot (float, m/s/day): Measured linear trend term
        gdot (float, m/s/day): Error on gdot
        gddot (float, m/s/day/day): Measured quadratic curvature 
                                    term
        gddot_err (float, m/s/day/day): Error on gddot
    
    Returns:
        log_likelihood_total (float): Likelihood of gdot AND gddot 
                                      conditioned on the model.
    """
    cdef double E, gdot_model, gddot_model
    cdef double log_likelihood_gdot, log_likelihood_gddot, log_likelihood_total
    
    E = M_2_evolvedE(M_anom_0, per, e, rv_epoch)

    gdot_model, gddot_model = gamma(a, m, e, i, om, E, per, m_star)

    # Log of the prefactor minus the log of the exponential
    log_likelihood_gdot  = log(1/(sqrt(two_pi)*gdot_err))\
                         - (gdot-gdot_model)**2/(2*gdot_err**2)
                         
    log_likelihood_gddot = log(1/(sqrt(two_pi)*gddot_err))\
                         - (gddot-gddot_model)**2/(2*gddot_err**2)
                 
    log_likelihood_total = log_likelihood_gdot + log_likelihood_gddot

    return log_likelihood_total


@cython.cdivision(True)
cdef (double, double) gamma(double a, double m, double e, 
                             double i, double om, double E, 
                             double per, double m_star):
    """
    Given an orbital model, calculates RV trend (gdot) and 
    curvature (gddot).

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
        gdot (float, m/s/day): Model linear trend term
        gddot (float, m/s/day/day): Model quadratic curvature term
    """

    cdef double sqrt_eterm,\
                sqrt_e_sq_term,\
                cos_E, sin_E,\
                nu, nu_dot, nu_ddot,\
                cos_nu_om, sin_nu_om, sin_i,\
                pre_fac, gamma_dot, gamma_ddot
             
    sqrt_eterm = sqrt((1+e)/(1-e))
    sqrt_e_sq_term = sqrt(1-e*e)
    
    cos_E = cos(E)
    sin_E = sin(E)

    #####################################
    ## Regular tan and arctan functions.
    ## Manipulations are possible with trig identities,
    ## but they don't give much of a speedup.
    nu = 2*atan(sqrt_eterm*tan(E/2))

    # nu derivatives use days (not seconds) to give gdot/gddot correct units 
    nu_dot = two_pi*sqrt_e_sq_term/(per*(1-e*cos_E)**2) # Units of day^-1
    nu_ddot = -nu_dot**2 * 2*e*sin_E/sqrt_e_sq_term # # Units of day^-2
    
    
    cos_nu_om = cos(nu+om)
    sin_nu_om = sin(nu+om)
    sin_i = sin(i)

    # Lovis+Fischer 2010 (analytic)
    pre_fac = sqrt(G)/sqrt_e_sq_term * m*sin_i/sqrt((m+m_star)*(a)) * auday2ms # AU/day ==> m/s

    gamma_dot = -pre_fac*nu_dot*sin_nu_om # m/s/day
    gamma_ddot = -pre_fac*(nu_dot**2*cos_nu_om + nu_ddot*sin_nu_om) # m/s/day/day


    return gamma_dot, gamma_ddot

    
cpdef M_2_evolvedE(double M0, double per, double e, double rv_epoch):
    """
    Use Radvel Kepler solver to take a mean anomaly at the beginning 
    of the Hip mission, evolve it to the RV epoch, and convert it to
    eccentric anomaly.
    
    Arguments:
        M0 (float, radians): Mean anomaly at the beginning of the 
                             Hipparcos mission
        per (float, days): Orbital period
        e (float): Orbital eccentricity
        rv_epoch (float, days): BJD at which the model RV trend and 
                                curv are calculated. This time should 
                                be close to the midpoint of the RV
                                baseline, where the measured trend
                                and curvature (presumably) fit the
                                time series best.
    
    Returns:
        E (float, radians): Eccentric anomaly at time = rv_epoch
  
    """
    # Fewer python references with no declarations. Not sure why.
    
    M_evolved = ((two_pi/per)*(rv_epoch - hip_beginning) + M0)%two_pi

    E = kepler_single(M_evolved, e)

    return E
    