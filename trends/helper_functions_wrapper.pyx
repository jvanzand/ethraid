# cython: language_level=3, boundscheck=False, cdivision=True, wraparound=False
## cython: binding=True
# The binding=True directive is useful for profiling code using kernprof. It also creates a lot of 'unused variable' warnings, so I will keep it commented out with two hashes for now
import os
import sys

path = os.getcwd()
sys.path.append(path+'/trends')

from kern_profiler_dummy import *

import numpy as np
cimport numpy as np
import scipy as sp
import scipy.stats as spst
from astropy.time import Time
import h5py

import cython
from libc.math cimport sin, cos, tan, atan, sqrt
from cpython cimport array

import radvel as rv
from c_kepler import _kepler as ck

##########################################
#### Kepler solver for one M and one e
# Wrapping kepler(M,e) a simple function that takes two doubles as
# arguments and returns a double
cdef extern from "../c_kepler/kepler.c":
    double kepler(double M, double e)
    double rv_drive(double t, double per, double tp, double e, double cosom, double sinom, double k )
##########################################
## Constants ##

cdef float pi, two_pi, math_e, G, M_sun, M_jup, au
cdef float hip_times[2]
cdef float gaia_times[2]

pi = 3.141592653589793
two_pi = 6.283185307179586
math_e = 2.718281828459045
G = 6.674299999999999e-08
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
au = 14959787070000.0

hip_times  = [Time(1989.85, format='decimalyear').jd, 
              Time(1993.21, format='decimalyear').jd] #https://www.cosmos.esa.int/web/hipparcos/catalogue-summary

gaia_times = [Time('2014-07-25', format='isot').jd, 
              Time('2017-05-28', format='isot').jd] #https://www.cosmos.esa.int/web/gaia/earlydr3

#@profile
def make_arrays(double m_star, tuple a_lim, tuple m_lim, double rv_epoch, int grid_num, int num_points):

    cdef double tp, a_min, a_max, m_min, m_max


    cdef np.ndarray[double, ndim=1] a_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    m_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    per_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    e_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    cosi_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    i_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    M_anom_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    T_anom_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    om_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    a_bins = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    m_bins = np.ndarray(shape=(num_points,), dtype=np.float64)

    cdef long [:] a_inds, m_inds

    np.random.seed(0)
    tp = 0
    a_min = a_lim[0]
    a_max = a_lim[1]
    m_min = m_lim[0]
    m_max = m_lim[1]


    # These semimajor axes are distances between the planet and the barycenter of the system. The star is on its own orbit, which we will get later.
    a_list = spst.loguniform.rvs(a_min, a_max, size=num_points)
    m_list = spst.loguniform.rvs(m_min, m_max, size=num_points)

    # Match up a_list and m_list and get the period for each pair (in days).
    per_list = P(a_list, m_star+m_list*(M_jup/M_sun) ) # Use this line when we are actually sampling a_tot, not a_planet
    #per_list = P(a_list * (m_star+m_list*(M_jup/M_sun))/m_star, m_star+m_list*(M_jup/M_sun) )

    # Eccentricities drawn from a beta distribution. I am using (a,b) = (0.867, 3.03) according to Winn & Fabrycky (2014).
    e_list = spst.beta(0.867, 3.03).rvs(num_points)
    e_list = np.where(e_list > 0.99, 0.99, e_list) # Replace e > 0.99 with 0.99

    cosi_list = np.random.uniform(0, 1, num_points)
    i_list = np.arccos(cosi_list)


    # Mean anomaly, uniformly distributed. This represents M at the beginning of the Hipparcos epoch for BOTH RVs and astrometry. Use this to solve for True anomaly.
    M_anom_0 = np.random.uniform(0, two_pi, num_points)

    # Evolving M_anom forward to the epoch of RV calculations.
    M_anom_evolved = M_anom_0 + two_pi*((rv_epoch - hip_times[0])/per_list)
    E_anom_rv = ck.kepler_array(M_anom_evolved, e_list) # Used in post_rv

    # Not evolving for use in astrometry calculations (b/c these are the STARTING angles ~1991)
    #E_anom_astro = ck.kepler_array(M_anom_list, e_list)
    #T_anom_astro = 2*np.arctan(np.sqrt((1+e_list)/(1-e_list)) * np.tan(E_anom_astro/2)) # Used in post_astro as T_anom_0.

    # Arguments of peri, uniformly distributed
    om_list = np.random.uniform(0, two_pi, num_points)

    # Longitudes of ascending node, uniformly distributed
    # Om_list = np.random.uniform(0, 2*pi, num_points)

    # Breaking up the (a, M) parameter space into grid_num x grid_num
    a_bins = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
    m_bins = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)

    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)

    return a_list, m_list, per_list, e_list, i_list, om_list, M_anom_0, E_anom_rv, a_inds, m_inds



cdef P(double [:] a, double [:] Mtotal):
    """
    Uses Kepler's third law to find the period of a planet (in days) given its
    semimajor axis and the total mass of the system.

    a (au): semi-major axis
    Mtotal (Msun): Mass of star + mass of object
    """
    cdef int size, i
    cdef double sec_2_days

    size = a.shape[0]
    sec_2_days = 1./(24*3600) # Note the 1.; with 1, the result would be 0 for some reason

    cdef np.ndarray[double, ndim=1] P_days = np.ndarray(shape=(size,), dtype=np.float64)

    for i in range(size):

        P_days[i] = sqrt((two_pi)**2*(a[i]*au)**3/(G*(Mtotal[i]*M_sun))) * sec_2_days

    return P_days

#@profile
def gamma_array(double m_star, double [:] a, double [:] Mp,
          double [:] per, double [:] e,
          double [:] i, double [:] om,
          double [:] E_anom):
    """
    Outsources intensive calculations to the pure-cython gamma function.
    """
    cdef int size, j

    size = a.shape[0]


    cdef np.ndarray[double, ndim=1]  gamma_dot = np.ndarray(shape=(size,), dtype=np.float64),\
                                    gamma_ddot = np.ndarray(shape=(size,), dtype=np.float64)


                                    
    for j in range(size):
       gamma_dot[j], gamma_ddot[j]  = gamma(m_star, a[j], Mp[j], per[j], e[j], i[j], om[j], E_anom[j])


    return gamma_dot, gamma_ddot

# It seems gamma() needs to be a cdef function, otherwise it returns nans
#@profile
cdef (double, double) gamma(double m_star, double a, double Mp, double per, 
      double e, double i, double om, double E):

    cdef double     a_cm, sqrt_eterm, tan_E2, nu,\
                    cos_E, tan_nu2, cos_E2, sin_i, cos_nu,\
                    sin_nu, cos_nu_om, sin_nu_om, sin_E,\
                    E_dot, nu_dot, prefac, gd_t1, gd_t2,\
                    gamma_dot, gamma_ddot

    Mp_g = Mp*M_jup
    m_star_g = m_star*M_sun
    
    #############
    a_cm = a*au
    #############
    
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
    pre_fac = sqrt(G)/sqrt_e_sq_term * Mp_g*sin_i/sqrt((Mp_g+m_star_g)*(a_cm)) * 1/100 # cm/s ==> m/s
    # Fischer (calculated, a)
    #pre_fac = 28.4329/sqrt_e_sq_term * Mp*sin_i/sqrt(m_star+Mp_g/M_sun)*1/sqrt(a_tot)
    # Fischer (calculated, per)
    #pre_fac = 28.4329/sqrt_e_sq_term * Mp*sin_i/(pow(m_star+Mp_g/M_sun, 0.6666666666)*pow(per/365.25, 0.3333333333))
    #print(pow(per/365.25, 0.333333333))
    
    gamma_dot = -pre_fac*nu_dot*sin_nu_om # m/s/d
    gamma_ddot = -pre_fac*(nu_dot**2*cos_nu_om + nu_ddot*sin_nu_om) # m/s/d/d

    return gamma_dot, gamma_ddot


#@profile
def rv_post(double gammadot, double gammadot_err,
                       double gammaddot, double gammaddot_err,
                       double m_star, double [:] a_list, double [:] m_list,
                       double [:] per_list, double [:] e_list, double [:] i_list,
                       double [:] om_list, double [:] E_anom_list, int num_points,
                       int grid_num):

    cdef double [:] m_tot_list = np.zeros(shape=(num_points), dtype=np.float64)
    cdef double [:] gammadot_list = np.zeros(shape=(num_points), dtype=np.float64)
    cdef double [:] gammaddot_list = np.zeros(shape=(num_points), dtype=np.float64)


    cdef int i, size
    cdef double chi_sq

    cdef double [:] rv_prob_list = np.zeros(shape=(num_points,), dtype=np.float64)
    #cdef double [:,:] rv_prob_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)
    cdef double model[2]
    cdef double data[2]
    cdef double err[2]


    gammadot_list, gammaddot_list = gamma_array(m_star, a_list, m_list, 
                                                per_list, e_list, i_list, om_list, E_anom_list)


    data[0] = gammadot
    data[1] = gammaddot

    err[0] = gammadot_err
    err[1] = gammaddot_err

    for i in range(num_points):

        model[0] = gammadot_list[i]
        model[1] = gammaddot_list[i]

        prob = likelihood_rv(model, data, err)

        rv_prob_list[i] = prob

    return rv_prob_list


#@profile
def astro_post(double delta_mu, double delta_mu_err, double m_star, double d_star,
               np.ndarray[double, ndim=1] a_list, double [:] m_list, double [:] per_list,
               double [:] e_list, double [:] i_list, double [:] om_list, double [:] M_anom_0_list,
               int num_points, int grid_num):


    #cdef double rot_mtrx_list[3][3]
    cdef double [:,::1] rot_mtrx # This makes rot_mtrx a memview
    rot_mtrx = np.zeros((3,3),dtype=np.float64)

    cdef double M_anom, E_anom, T_anom, r_pl, r_star, mass_ratio

    cdef double r_vec_list[3]
    cdef double [:] r_vec = r_vec_list

    cdef double rot_vec_list[3]
    cdef double [:] rot_vec = rot_vec_list

    cdef double ang_pos_list[2]
    cdef double [:] ang_pos = ang_pos_list

    ###############################################
    cdef int j, k, l, a_j, m_j

    cdef double time_steps[2]
    cdef double both[2][2]

    cdef double baseline_yrs, start_time, end_time, elapsed_time, a, m, per, e, i, om, M_anom_0
    cdef double mass_ratio_constant, cm_2_mas, cms_2_masyr
    cdef double mean_motion, sqrt_eterm, a_units, e_sq, r_star_num_fac
    cdef double M1, M2, E1, E2, x_pos_avg, y_pos_avg, x_vel_avg, y_vel_avg
    cdef double delta_mu_model, chi_sq, prob

    cdef double time_endpoints[2][2]
    cdef double ang_pos_avg[2][2]
    cdef double mu_avg[2][2]
    cdef double mu_gaia[2]
    cdef double mu_hg[2]

    # I don't think I can declare a c-type array with size (grid_num, grid_num)
    cdef double [:] astro_prob_list = np.zeros(shape=(num_points), dtype=np.float64)
    cdef double [:,:] astro_prob_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)

    # Time in days between epoch mid-points
    baseline_yrs = ((gaia_times[1] + gaia_times[0])/2 - (hip_times[1] + hip_times[0])/2)/365.25

    cm_2_mas = (206265*1e3)/d_star
    cmd_2_masyr = cm_2_mas * 365.25 # cm/day to milli-arcseconds/year

    time_endpoints = [[hip_times[0], gaia_times[0]], [hip_times[1], gaia_times[1]]]

    ##############################################

    cdef double v_vec_pl_list[3]
    cdef double [:] v_vec_pl = v_vec_pl_list

    cdef double v_vec_star_list[3]
    cdef double [:] v_vec_star = v_vec_star_list

    cdef double rotated_v_vec_list[3]
    cdef double [:] rotated_v_vec = rotated_v_vec_list

    cdef double mu[2]
    ###############################################

    for j in range(num_points): # Loop over the desired number of random points

        a = a_list[j]
        m = m_list[j]
        per = per_list[j]
        e = e_list[j]
        i = i_list[j]
        om = om_list[j]
        M_anom_0 = M_anom_0_list[j]
        
        a_units = a*au
        #a_star_units = a_units*mass_ratio # Calculate the sma (cm) of the star's orbit around barycenter. 
        
        ######################
        mass_ratio = (m*M_jup/(m_star*M_sun + m*M_jup))
        a_star_units = a_units*mass_ratio
        ######################
        
        mean_motion = two_pi/per # Units of 1/day because per is in days
        sqrt_eterm = sqrt((1+e)/(1-e))
        e_sq = e**2
        rot_matrix(i, om, 0, rot_mtrx) # Omega = 0 arbitrarily

        for l in range(2): # Hipparcos or Gaia
            start_time = time_endpoints[0][l] - time_endpoints[0][0] # The "start time" of Hip or Gaia relative to the start of Hip. For Hip, start_time is 0. For Gaia, it is the time between Hip_start and Gaia_start
            end_time = time_endpoints[1][l] - time_endpoints[0][0] # End time relative to the start of Hip.
            
            
            ## Mean anomaly is the elapsed time times the mean motion, plus a randomly-sampled starting mean anomaly
            M1 = mean_motion*start_time + M_anom_0
            M2 = mean_motion*end_time + M_anom_0
            
            E1 = kepler(M1, e)
            E2 = kepler(M2, e)
            
            # Get position of the STAR.
            x_pos_avg, y_pos_avg = pos_avg(a_star_units, mean_motion, e, E1, E2, start_time, end_time)
            

            # r_vec points from barycenter to the *star* (note the - sign) in the orbital plane, and has magnitude r_star. Like r_star, it has units of cm.
            # Since we're using the E_anom for the planet, the star is located in the opposite direction
            r_vec[0] = -x_pos_avg
            r_vec[1] = -y_pos_avg
            r_vec[2] = 0
            

            # rot_vec points from barycenter to the star, but in coordinates where the xy-plane is the sky plane and the z-axis points toward Earth
            mat_mul(rot_mtrx, r_vec, rot_vec)


            # ang_pos is the angular position of the star relative to barycenter in milli-arcseconds.
            ang_pos[0] = rot_vec[0]*cm_2_mas
            ang_pos[1] = rot_vec[1]*cm_2_mas
            
            
            ang_pos_avg[l][0] = ang_pos[0]
            ang_pos_avg[l][1] = ang_pos[1]
            
            ################### Angular Velocities ########################
            
            # I only need Gaia velocities, not Hip velocities
            if l == 0:
                continue
            
            # Get velocity of the star
            x_vel_avg, y_vel_avg = vel_avg(a_star_units, mean_motion, e, E1, E2, start_time, end_time)

            # Since we're using the E_anom for the planet, the star is moving in the opposite direction
            v_vec_star[0] = -x_vel_avg
            v_vec_star[1] = -y_vel_avg
            v_vec_star[2] = 0


            mat_mul(rot_mtrx, v_vec_star, rotated_v_vec)

            # mu is the proper motion of the star due to the planet's orbit in milli-arcseconds per year.
            # Since period is in days and a_star_units in cm, velocities are in cm/day.
            mu[0] = rotated_v_vec[0]*cmd_2_masyr
            mu[1] = rotated_v_vec[1]*cmd_2_masyr
            
            # mu_avg is a 2x2 array. The top row stays empty because we skip Hip. The bottom row is Gaia prop. motion
            mu_avg[l][0] = mu[0]
            mu_avg[l][1] = mu[1]

            ###############################################################
            ###############################################################


        mu_gaia = mu_avg[1]

        # To get the positional avg., subtract the epoch positions and divide by the time between in years.
        # First index tells Hip ([0]) or Gaia ([1]), second index tells x ([0]) or y ([1])
        # Units of mas/yr
        mu_hg[0] = (ang_pos_avg[1][0] - ang_pos_avg[0][0])/baseline_yrs # x-comp. = gaia_x - hip_x
        mu_hg[1] = (ang_pos_avg[1][1] - ang_pos_avg[0][1])/baseline_yrs # y-comp. = gaia_y - hip_y


        delta_mu_model = sqrt((mu_hg[0] - mu_gaia[0])**2 + (mu_hg[1] - mu_gaia[1])**2)

        prob = likelihood_astro(delta_mu, delta_mu_err, delta_mu_model)

        astro_prob_list[j] = prob

    return astro_prob_list

cdef pos_avg(double a, double n, double e, double E1, double E2, 
                  double t1, double t2):
    """
    Calculate the average x/y positions of an object on an elliptical orbit, 
    where (0,0) is the focus.
    
    a (cm): semi-major axis
    n (1/days): 2pi/per
    t1, t2 (days): beginning and ending time to calculate average
    
    returns: Average x and y positions (cm)
    """
    
    cdef double x_term_1, x_term_2, x_integral, x_avg,\
                y_term_1, y_term_2, y_integral, y_avg
    
    x_term_1 = a/n * ((1+e**2)*sin(E1) - e/4 * (6*E1+sin(2*E1)))
    x_term_2 = a/n * ((1+e**2)*sin(E2) - e/4 * (6*E2+sin(2*E2)))
    
    x_integral = x_term_2 - x_term_1
    
    x_avg = 1/(t2-t1) * x_integral
    
    y_term_1 = -a*sqrt(1-e**2)/n * (cos(E1) - e/2 * cos(E1)**2 )
    y_term_2 = -a*sqrt(1-e**2)/n * (cos(E2) - e/2 * cos(E2)**2 )
     
    y_integral = y_term_2 - y_term_1
     
    y_avg = 1/(t2-t1) * y_integral
    
    
    return x_avg, y_avg

cdef vel_avg(double a, double n, double e, double E1, double E2, 
                  double t1, double t2):
    """
    Calculate the average x/y positions of an object on an elliptical orbit, 
    where (0,0) is the focus.

    a (cm): semi-major axis
    n (1/days): 2pi/per
    !!! THHIS FUNCTION doesn't actually use n. It can be removed.
    t1, t2 (days): beginning and ending time to calculate average
    
    returns: Average x and y velocities (cm/day)
    """
    
    cdef double x_term_1, x_term_2, x_integral, x_avg,\
                y_term_1, y_term_2, y_integral, y_avg,\
                r1, r2
    
    r1 = a*(1-e*cos(E1))
    r2 = a*(1-e*cos(E2))

    x_term_1 = a**2/r1 * (cos(E1) + e/2 * sin(E1)**2)
    x_term_2 = a**2/r2 * (cos(E2) + e/2 * sin(E2)**2)

    x_integral = x_term_2 - x_term_1

    x_avg = 1/(t2-t1) * x_integral


    y_term_1 = a**2*sqrt(1-e**2)/r1 * (sin(E1) - e/2 * (E1 + 0.5*sin(2*E1)) )
    y_term_2 = a**2*sqrt(1-e**2)/r2 * (sin(E2) - e/2 * (E2 + 0.5*sin(2*E2)) )
 
    y_integral = y_term_2 - y_term_1
 
    y_avg = 1/(t2-t1) * y_integral


    return x_avg, y_avg


cdef double likelihood_astro(double data, double data_err, double model):
    """
    Simple calculation of likelihood given data and model
    assuming gaussian uncertainties and uniform priors.
    """

    cdef double chi_sq, prob


    chi_sq = ((data-model)/data_err)**2

    prob = math_e**(-chi_sq/2)

    return prob

cdef double likelihood_rv(double [:] model, double [:] data, double [:] data_err):
    """
    Simple calculation of likelihood given data and model
    assuming gaussian uncertainties and uniform priors. Accepts lists of
    data points because RV uses both gdot and gddot.
    """

    cdef int i, size
    cdef double chi_sq, prob

    size = data.shape[0]
    chi_sq = 0.0

    for i in range(size):

        chi_sq += ((data[i]-model[i])/data_err[i])**2

    prob = math_e**(-chi_sq/2)

    return prob

def prob_array(double [:] prob_list, long [:] a_inds, long [:] m_inds, int grid_num):
    """
    Form a list of probabilities into a 2D array
    """

    cdef int i, size, a_i, m_i
    cdef double [:,:] prob_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)

    size = prob_list.shape[0]

    for i in range(size):

        a_i = a_inds[i]
        m_i = m_inds[i]

        prob_array[m_i, a_i] += prob_list[i]

    return np.array(prob_array)



def post_tot(double [:] rv_post_list, double [:] astro_post_list, int grid_num,
            long [:] a_inds, long [:] m_inds):

    cdef int size, i, a_i, m_i
    cdef double [:,:] tot_prob_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)
    cdef double prob

    size = rv_post_list.size

    for i in range(size):

        a_i = a_inds[i]
        m_i = m_inds[i]

        prob = rv_post_list[i]*astro_post_list[i]

        tot_prob_array[m_i, a_i] += prob

    return np.array(tot_prob_array)

#@profile
cdef void rot_matrix(double i, double om, double Om, double [:,::1] rot_mtrx):
    """
    This is P3*P2*P1 from Murray & Dermott. It is not given explicitly in the text. 
    They multiply it immediately by r*[cos(f), sin(f), 0] because this gives the 
    projection of position onto the sky. However, we also need the projection of 
    velocity, so we need the matrix before multiplication by the position vector.

    This function doesn't return anything. Instead, declare a matrix in your 
    function and this will update it, saving lots of time by not allocating memory 
    to and returning a matrix.
    """
    cdef double sin_Om, sin_om, sin_i, cos_Om, cos_om, cos_i

    sin_Om = sin(Om)
    sin_om = sin(om)
    sin_i  = sin(i)
    cos_Om = cos(Om)
    cos_om = cos(om)
    cos_i  = cos(i)
    sin_Om_cos_i = sin_Om*cos_i
    cos_Om_cos_i = cos_Om*cos_i

    rot_mtrx[0][0] = cos_Om*cos_om - sin_Om_cos_i*sin_om
    rot_mtrx[0][1] = -sin_om*cos_Om - sin_Om_cos_i*cos_om
    rot_mtrx[0][2] = sin_Om*sin_i

    rot_mtrx[1][0] = sin_Om*cos_om + cos_Om_cos_i*sin_om
    rot_mtrx[1][1] = -sin_om*sin_Om + cos_Om_cos_i*cos_om
    rot_mtrx[1][2] = -cos_Om*sin_i

    rot_mtrx[2][0] = sin_i*sin_om
    rot_mtrx[2][1] = sin_i*cos_om
    rot_mtrx[2][2] = cos_i

    #return rot_mtrx


#@profile
cdef void mat_mul(double [:,:] mat, double [:] in_vec, double [:] out_vec):
    """
    This is written specifically to matrix multiply rot_matrix (3x3) with
    r_unit_vec (3x1) and later v_vec_star (3x1) in astro_post_dense_loop.
    """

    cdef int i, k

    for i in range(3):
        out_vec[i] = 0
        for k in range(3):
            out_vec[i] += mat[i][k]*in_vec[k]


def contour_levels(prob_array, sig_list, t_num = 1e3):
    """
    Contour drawing method taken from https://stackoverflow.com/questions/37890550/python-plotting-percentile-contour-lines-of-a-probability-distribution
    This function takes a 2-D array of probabilities and returns a 1-D array of the probability values corresponding to 1-sigma and 2-sigma
    contours. In this case, the 1-sigma contour encloses 68% of the total probability. The array is expected to be normalized. sig_list is
    a list containing any combination of the integers 1, 2, or 3 to indicate desired contours. For example, [1,3] will return the 1 and 3
    sigma contours.
    This function uses scipy.interpolate.interp1d.
    """


    # An array of probabilites from 0 to prob_max in rate_array
    t = np.linspace(0, np.array(prob_array).max(), int(t_num))

    # (prob_array >= t[:, None, None]) is a 3D array of shape (array_num, array_num, t_num). Each (array_num, array_num) layer is a 2D array of bool values indicating which values are greater than the value of the given t step.
    # Multiplying this 3D array of bools by prob_array replaces the bools with the array value if the bool is T and 0 if the bool is F.
    # Finally, sum along the array_num axes to get a single list of values, each with the total summed probability in its array.

    # integral is a 1D array of floats. The ith float is the sum of all probabilities in prob_array greater than the ith probability in t

    integral = ((prob_array >= t[:, None, None])*prob_array).sum(axis=(1,2))

    # Now create a function that takes integral as the x (not the y) and then returns the corresponding prob value from the t array. Interpolating between integral values allows me to choose any enclosed total prob. value (ie, integral value) and get the corresponding prob. value to use as my contour.
    f = sp.interpolate.interp1d(integral, t)

    contour_list = []
    prob_list = [0.68, 0.95, 0.997]

    for i in sig_list:
        contour_list.append(prob_list[i-1])

    # The plt.contourf function requires at least 2 levels. So if we want just one level, include a tiny contour that encompasses a small fraction of the total probability.
    if len(sig_list) == 1:
        contour_list.append(contour_list[0]-1e-4)
        # contour_list.append(1e-3)

    # Make sure list is in descending order
    t_contours = f(np.array(sorted(contour_list, reverse=True)))

    return t_contours

def contour_levels_1D(prob_list, sig_list, t_num = 1e3):
    """
    Same as contour_levels, but adapted for 1D arrays. Hopefully I can condense these into 1 in the future.
    """


    # An array of probabilites from 0 to prob_max in rate_array
    t = np.linspace(0, prob_list.max(), int(t_num))

    # integral is a 1D array of floats. The ith float is the sum of all probabilities in prob_array greater than the ith probability in t

    integral = ((prob_list >= t[:, None])*prob_list).sum(axis=(1))

    # Now create a function that takes integral as the x (not the y) and then returns the corresponding prob value from the t array. Interpolating between integral values allows me to choose any enclosed total prob. value (ie, integral value) and get the corresponding prob. value to use as my contour.
    f = sp.interpolate.interp1d(integral, t)

    contour_list = []
    prob_list = [0.68, 0.95, 0.997]

    for i in sig_list:
        contour_list.append(prob_list[i-1])

    # The plt.contourf function requires at least 2 levels. So if we want just one level, include a tiny contour that encompasses a small fraction of the total probability. In this case, the contour we want will be at the 0th index.
    if len(sig_list) == 1:
        contour_list.append(contour_list[0]-1e-4)

    # Make sure the list of integrals is in descending order (eg, 99.7%, 95%, 68%). This will make the list of probabilities be in ascending order (eg, 0.05, 0.01, 0.007). These correspond to descending sigma levels (3, 2, 1).
    t_contours = f(np.array(sorted(contour_list, reverse=True)))


    return t_contours

def bounds_1D(prob_array, value_spaces, interp_num = 1e4):
    """
    Given a 2D probability array, this function collapses the array along each 
    axis to find the 68% confidence interval.

    value_spaces represents the parameter intervals covered by the array along each axis.
    It is expected in the form [(min_value1, max_value1), (min_value2, max_value2)], 
    where 1 and 2 refer to the 0th and 1st axes.
    Note that the limits MUST be in this order: if the array has shape (x_num, y_num), 
    then value_spaces must be [x_lims, y_lims].
    """
    lvls_2sig_list = []
    inds_2sig_list = []
    bounds_list = []
    
    # First compute bounds for a by collapsing along the m (aka i=0) axis
    for i in range(2):

        array_1D = prob_array.sum(axis=i)
        grid_num = len(array_1D)


        # This gives only the 2-sigma contour level, so that we get the 2-sigma limits at the end
        lvl_2sig = contour_levels_1D(array_1D, [2])[0]

        # Interpolate between the points to get a finer spacing of points. This allows for more precise parameter estimation.
        func = sp.interpolate.interp1d(range(grid_num), array_1D)

        # Array over the same interval, but spaced (probably) more finely
        fine_array = np.linspace(0, grid_num-1, int(interp_num))

        # This is analogous to the original array_1D, but finer
        interp_vals = func(fine_array)

        #import matplotlib.pyplot as plt

        #plt.plot(range(len(fine_array)), interp_vals)
        #plt.show()


        # This is a shaky step. I'm just looking for places where the function value is really close to the probability corresponding to 2-sigma. But from what I can tell, this will fall apart for multimodal distributions, and maybe in other cases too. I use the 'take' method to pick out the first and last indices.


        inds_2sig = np.where(abs(interp_vals - lvl_2sig) < 1e-2*lvl_2sig)[0].take((0,-1))

        # value_bounds is a tuple of actual values, not indices
        value_bounds = index2value(inds_2sig, (0, interp_num-1), value_spaces[::-1][i])
        
        lvls_2sig_list.append(lvl_2sig)
        inds_2sig_list.append(inds_2sig * grid_num/interp_num)
        
        bounds_list.append(value_bounds)

    return bounds_list, lvls_2sig_list, inds_2sig_list


def value2index(value, index_space, value_space):
    """
    The inverse of index2value: take a value on a
    log scale and convert it to an index. index_space
    and value_space are expected as tuples of the form
    (min_value, max_value).
    """

    min_index, max_index = index_space[0],  index_space[1]
    min_value, max_value = value_space[0], value_space[1]
    
    value = np.array(value)
    value = value[min_value <= value <= max_value+0.1] # Add 0.1 because highest value was getting cut off

    index_range = max_index - min_index
    log_value_range = np.log10(max_value) - np.log10(min_value)
    

    index = (np.log10(value)-np.log10(min_value))*(index_range/log_value_range) + min_index
    
    # Looking back at this, I don't see any reason to round these values to integers. I am converting them into "index space," but they don't need to be integers for me to plot them.
    #if isinstance(value, (int, float)):
    #    ind = int(np.around(index))
    #else:
    #    ind = [int(np.around(i)) for i in index]

    return index

def index2value(index, index_space, value_space):
    """
    The axis values for a plotted array are just the array indices.
    I want to convert these to Msini and a values, and on a log
    scale. This function takes a single index from a linear index range,
    and converts it to a parameter value in log space. index_space and
    value_space are expected as tuples of the form (min_value, max_value).
    index is in the range of index_space.
    """
    index = np.array(index)

    min_index, max_index = index_space
    min_value, max_value = value_space

    index_range = max_index - min_index
    log_value_range = np.log10(max_value) - np.log10(min_value)

    # Convert from a linear space of indices to a linear space of log(values).
    log_value = (index-min_index)*(log_value_range/index_range) + np.log10(min_value)

    value = np.around(10**(log_value), 2) # Round to 2 decimal places

    return value

def period_lines(m, per, m_star):
    """
    Function to draw lines of constant period on the final plot.
    Rearranges Kepler's 3rd law to find how semi-major axis a 
    varies with period, companion mass, and stellar mass.
    
    Intended usage: Calculate an array of a values for a fixed per
                and m_star and an array of companion masses.
                
    Arguments:
            m (list of floats): companion masses (M_J)
            per (float): companion orbital period (days)
            m_star (float): stellar mass (M_sun)
    
    Returns:
            a (list of floats): Semi-major axis values (au)
    """
    m_grams = m*M_jup
    per_sec = per*24*3600
    m_star_grams = m_star*M_sun
    
    a_cm = ((per_sec/two_pi)**2*G*(m_grams+m_star_grams))**(0.3333333333)
    
    
    return a_cm / au #* (m_star_grams/(m_grams+m_star_grams))
    
    
def save_array(array, filename):

    h5f = h5py.File(filename, 'w')
    h5f.create_dataset('dataset_1', data=array)
    h5f.close()
    
    return









