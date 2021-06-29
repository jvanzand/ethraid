# cython: language_level=3, boundscheck=False, cdivision=True, wraparound=False
from kern_profiler_dummy import *

import numpy as np
cimport numpy as np
import scipy.stats as spst
import cython
from libc.math cimport sin, cos, tan, atan, sqrt
from cpython cimport array

#from scipy.stats import loguniform, beta

import scipy as sp
#import astropy.constants as c
#import astropy.units as u

# Kind of sloppy. This is just so you can access c_kepler from no_classes/
import sys
sys.path.append('../')

import radvel as rv

from c_kepler import _kepler as ck
import helper_functions as hlp

##########################################
#### Kepler solver for one M and one e 
# Wrapping kepler(M,e) a simple function that takes two doubles as
# arguments and returns a double
cdef extern from "../c_kepler/kepler.c":
    double kepler(double M, double e)
    double rv_drive(double t, double per, double tp, double e, double cosom, double sinom, double k )
##########################################
## Constants ##

cdef float pi, G, M_sun, M_jup, au

pi = 3.141592653589793
math_e  = 2.718281828459045
G =  6.674299999999999e-08
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
au = 14959787070000.0


def make_arrays(double m_star, tuple a_lim, tuple m_lim, int grid_num, int num_points):

    cdef double tp, a_min, a_max, m_min, m_max, two_pi

    
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
                                
    cdef long [:]                   a_inds, m_inds
    
    np.random.seed(0)
    tp = 0
    a_min = a_lim[0]
    a_max = a_lim[1]
    m_min = m_lim[0]
    m_max = m_lim[1]
    two_pi = 2*pi

    
    # These semimajor axes are distances between the planet and the barycenter of the system. The star is on its own orbit, which we will get later.
    a_list = spst.loguniform.rvs(a_min, a_max, size=num_points)
    m_list = spst.loguniform.rvs(m_min, m_max, size=num_points)
    
    # Match up a_list and m_list and get the period for each pair (in days).
    per_list = P(a_list, m_star+m_list*(M_jup/M_sun) )
    
    # Eccentricities drawn from a beta distribution. I am using (a,b) = (0.867, 3.03) according to Winn & Fabrycky (2014)
    e_list = spst.beta(0.867, 3.03).rvs(num_points)
    

    cosi_list = np.random.uniform(0, 1, num_points)
    i_list = np.arccos(cosi_list)
    #sini_list = np.sqrt(1-cosi_list**2)
    
    
    # Mean anomaly, uniformly distributed. Use this to solve for True anomaly.
    M_anom_list = np.random.uniform(0, two_pi, num_points)

    # Solving Kepler now to use in RV posterior
    E_anom_list = ck.kepler_array(M_anom_list, e_list)
    T_anom_list = 2*np.arctan(np.sqrt((1+e_list)/(1-e_list)) * np.tan(E_anom_list/2))

    # Arguments of peri, uniformly distributed
    om_list = np.random.uniform(0, two_pi, num_points)
    
    # Longitudes of ascending node, uniformly distributed
    #Om_list = np.random.uniform(0, 2*pi, num_points)

           
    # Breaking up the (a, M) parameter space into grid_num x grid_num
    a_bins = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
    m_bins = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)
    
    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)
    
    
    
    return a_list, m_list, per_list, e_list, i_list, om_list, M_anom_list, E_anom_list, T_anom_list, a_inds, m_inds



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
    sec_2_days = 1./(24*3600) # Note the 1.; without it the result would be 0

    cdef np.ndarray[double, ndim=1] P_days = np.ndarray(shape=(size,), dtype=np.float64)
    
    for i in range(size):
    
        P_days[i] = sqrt((2*pi)**2*(a[i]*au)**3/(G*(Mtotal[i]*M_sun))) * sec_2_days
    
    return P_days
    

def gamma_array(double [:] a, double [:] Mp, 
          double [:] per, double [:] e, 
          double [:] i, double [:] om, 
          double [:] E_anom):
    """
    Outsources intensive calculations to the pure-cython gamma function. 
    Unfortunately, this is slower than gamma_direct_FAST below, which uses numpy and runs in ~0.5x the time.
    I think I need to import the functions from helper_functions.c, but this throws errors that I haven't figured out.
    """
    print('Using gamma_array')
    cdef int size, j
    
    size = a.shape[0]
    
    
    cdef np.ndarray[double, ndim=1]  gamma_dot = np.ndarray(shape=(size,), dtype=np.float64),\
                                    gamma_ddot = np.ndarray(shape=(size,), dtype=np.float64)

    for j in range(size):
       gamma_dot[j], gamma_ddot[j]  = gamma(a[j], Mp[j], per[j], e[j], i[j], om[j], E_anom[j])
        

    return gamma_dot, gamma_ddot


cdef (double, double) gamma(double a, double Mp, double per, double e, double i, double om, double E):

    cdef double Mp_units, a_units, sqrt_eterm, tan_E2, nu,\
                cos_E, tan_nu2, cos_E2, sin_i, cos_nu,\
                sin_nu, cos_nu_om, sin_nu_om, sin_E,\
                E_dot, nu_dot, prefac, gd_t1, gd_t2,\
                gamma_dot, gd_t1_dot, gd_t2_dot, gdd_t1, gdd_t2, gamma_ddot

    Mp_units = Mp*M_jup
    a_units = a*au
    
    sqrt_eterm = sqrt((1+e)/(1-e))
    tan_E2 = tan(E/2)
    #nu = 2*atan(sqrt((1+e)/(1-e))**0.5*tan(E/2))
    nu = 2*atan(sqrt_eterm*tan_E2)
    
    
    cos_E = cos(E)
    tan_nu2 = tan(nu/2)
    cos_E2 = cos(E/2)
    cos_E2_2 = 1/(cos_E2*cos_E2)
    sin_i = sin(i)
    
    cos_nu = cos(nu)
    sin_nu = sin(nu)
    cos_nu_om = cos(nu+om)
    sin_nu_om = sin(nu+om)
    sin_E = sin(E)

    # Differentiate Kepler's equation in time to get E_dot
    # Note that E_dot has units of (1/per), where [per] is days. Therefore [gamma_ddot] = m/s/d^2
    E_dot = (2*pi/per)/(1-e*cos_E)
    #nu_dot = (1+tan(nu/2)**2)**-1 * ((1+e)/(1-e))**0.5 * cos(E/2)**-2 * E_dot
    nu_dot = (1+tan_nu2**2)**-1 * sqrt_eterm * cos_E2_2 * E_dot

    # Convert prefac units from cm/s^2 to m/s/day
    # Negative just depends on choice of reference direction. I am being consistent with radvel rv_drive function.
    prefac = -(Mp_units*G*sin_i)/(a_units**2*(1-e)) * 864 # Save calculation of 24*3600 / 100


    gd_t1 = (1+cos_nu)/(1+cos_E)
    gd_t2 = sin_nu_om/(1-e*cos_E)


    gamma_dot = prefac*gd_t1*gd_t2

    gd_t1_dot = ((1+cos_nu)*sin_E * E_dot - (1+cos_E)*sin_nu*nu_dot) / (1+cos_E)**2
    gd_t2_dot = ((1-e*cos_E)*cos_nu_om * nu_dot - sin_nu_om*e*sin_E*E_dot) / (1-e*cos_E)**2


    gdd_t1 = gd_t2 * gd_t1_dot
    gdd_t2 = gd_t1 * gd_t2_dot

    gamma_ddot = prefac*(gdd_t1+gdd_t2)

    return gamma_dot, gamma_ddot
    
    
#@profile
def rv_post(double gammadot, double gammadot_err, 
                       double gammaddot, double gammaddot_err,
                       double m_star, double [:] a_list, double [:] m_list,
                       double [:] per_list, double [:] e_list, double [:] i_list, 
                       double [:] om_list, double [:] E_anom_list, int num_points,
                       int grid_num, long [:] a_inds, long [:] m_inds):
    
    cdef double [:] m_tot_list = np.zeros(shape=(num_points), dtype=np.float64)
    cdef double [:] gammadot_list = np.zeros(shape=(num_points), dtype=np.float64)
    cdef double [:] gammaddot_list = np.zeros(shape=(num_points), dtype=np.float64)
    
    
    cdef int i, size, a_i, m_i  #Typing a_i and m_i slows it down? Double check.
    cdef double chi_sq
    
    cdef double [:] rv_prob_list = np.zeros(shape=(num_points,), dtype=np.float64)
    cdef double [:,:] rv_prob_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)
    
    
    gammadot_list, gammaddot_list = gamma_array(a_list, m_list, per_list, e_list, i_list, om_list, E_anom_list)
    
    for i in range(num_points):
    
        a_i = a_inds[i]
        m_i = m_inds[i]
    
        chi_sq = ((gammadot_list[i]-gammadot)/(gammadot_err))**2 + ((gammaddot_list[i]-gammaddot)/(gammaddot_err))**2
        prob = math_e**(-chi_sq/2)
        
        rv_prob_array[m_i, a_i] += prob
        rv_prob_list[i] = prob
    
    return rv_prob_array, rv_prob_list

#@profile
def astro_post_array(double delta_mu, double delta_mu_err, double m_star, double d_star, 
                        np.ndarray[double, ndim=1] a_list, double [:] m_list, double [:] per_list,
                         double [:] e_list, double [:] i_list, double [:] om_list, double [:] T_anom_0_list, 
                         int num_points, int grid_num, long [:] a_inds, long [:] m_inds, int t_num):
    """
    M_anom_prog is not randomly-sampled. It is a deterministic list, based on elapsed time and per (which is sampled).
    It also has 2 dimensions: t_num and per_num (one array of mean anomalies through t_num for each per_num)
    The array of semi-major axes (a) has to be a np array to support scalar multiplication.
    """
    cdef int j, k, l, a_j, m_j
    

    cdef double hip_times[2]
    cdef double gaia_times[2]
    #cdef double time_endpoints[2][2]
    cdef double time_steps[2]
    cdef double both[2][2]

                                                         
    cdef double baseline_yrs, start_time, end_time, elapsed_time, a, m, per, e, i, om, T_anom_0, cm_2_mas, cms_2_masyr
    cdef double ang_pos_x_sum, ang_pos_y_sum, mu_x_sum, mu_y_sum
    cdef double delta_mu_model, chi_sq, prob
    
    cdef double time_endpoints[2][2]
    cdef double ang_pos_avg[2][2]
    cdef double mu_avg[2][2]
    cdef double mu_gaia[2]
    cdef double mu_hg[2]

    # I don't think I can declare a c-type array with size (grid_num, grid_num)
    cdef double [:] astro_prob_list = np.zeros(shape=(num_points), dtype=np.float64)
    cdef double [:,:] astro_prob_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)
    
    hip_times  = [2447837.75, 2449065.15]
    gaia_times = [2456863.5, 2457531.5]
    
    # Time in days between epoch mid-points
    baseline_yrs = ((gaia_times[1] + gaia_times[0])/2 - (hip_times[1] + hip_times[0])/2)/365

    
    cm_2_mas = (206265*1e3)/d_star
    cmd_2_masyr = cm_2_mas * 365.25 # cm/day to milli-arcseconds/year
    
    time_endpoints = [[hip_times[0], gaia_times[0]], [hip_times[1], gaia_times[1]]]
    
    # Best to calculate these now so they aren't re-calculated many times in each loop
    time_steps[0] = (hip_times[1] - hip_times[0])/(t_num+1)
    time_steps[1] = (gaia_times[1] - gaia_times[0])/(t_num+1)
    
    for j in range(num_points): # Loop over the desired number of random points
    
        a = a_list[j]
        m = m_list[j]
        per = per_list[j]
        e = e_list[j]
        i = i_list[j]
        om = om_list[j]
        T_anom_0 = T_anom_0_list[j]
        
        for l in range(2): # Hipparcos or Gaia
            start_time = time_endpoints[0][l] - time_endpoints[0][0] # The "start time" of Hip or Gaia relative to the start of Hip. For Hip, start_time is 0. For Gaia, it is the time between Hip_start and Gaia_start
            
            time_step = time_steps[l]
            
            ang_pos_x_sum = 0
            ang_pos_y_sum = 0
            
            mu_x_sum = 0
            mu_y_sum = 0
            
            for k in range(t_num+1): # Start at 0, finish at t_num. This is a loop over the time steps of each mission

                elapsed_time = k*time_step + start_time

                both = astro_post(m_star, cm_2_mas, cmd_2_masyr, a, m, per, e, i, om, T_anom_0, num_points, grid_num, a_inds, m_inds, t_num, elapsed_time)
                
                ang_pos = both[0]
                mu = both[1]
                
                # Sum up each coordinate independently in anticipation of taking an average
                ang_pos_x_sum += ang_pos[0]
                ang_pos_y_sum += ang_pos[1]
            
                mu_x_sum += mu[0]
                mu_y_sum += mu[1]
                
            
            # Once we have the sums over the whole epoch, divide by number of points to get the avg.
            ang_pos_avg[l][0] = ang_pos_x_sum/(t_num+1)
            ang_pos_avg[l][1] = ang_pos_y_sum/(t_num+1)            
        
            mu_avg[l][0] = mu_x_sum/(t_num+1)
            mu_avg[l][1] = mu_y_sum/(t_num+1)
            
            #print('horses', ang_pos_avg)
            #print(mu_avg)
            #dfdf
        
        # The gaia mu is just the 1th component. mu_gaia has dimension (2) and units of mas/yr
        mu_gaia = mu_avg[1]

        # To get the positional avg., subtract the epoch positions and divide by the time between in years.
        # Units of mas/yr
        mu_hg[0] = (ang_pos_avg[1][0] - ang_pos_avg[0][0])/baseline_yrs # x-comp. = gaia_x - hip_x
        mu_hg[1] = (ang_pos_avg[1][1] - ang_pos_avg[0][1])/baseline_yrs # y-comp. = gaia_y - hip_y
        

        delta_mu_model = sqrt((mu_hg[0] - mu_gaia[0])**2 + (mu_hg[1] - mu_gaia[1])**2)
            
        
        chi_sq = ((delta_mu-delta_mu_model)/delta_mu_err)**2
        
        
        #if j < 10:
        #    print('Yeeep')
        #    #print(delta_mu_model)
        #    print(chi_sq)
        #    print(-chi_sq/2)
        #    x = -chi_sq/2
        #    print(math_e**x)
        #    print(math_e**(-chi_sq/2))
        #    print('\n')
        #else:
        #    asdfasdf
        
        a_j = a_inds[j]
        m_j = m_inds[j]
        
        prob = math_e**(-chi_sq/2)
        
        astro_prob_list[j] = prob
        astro_prob_array[m_j, a_j] += prob
        
    
    return astro_prob_array, astro_prob_list

#@profile
def astro_post(double m_star, double cm_2_mas, double cmd_2_masyr, double a, double m, double per, 
               double e, double i, double om, double T_anom_0, int num_points, int grid_num, 
               long [:] a_inds, long [:] m_inds, int t_num, double elapsed_time):

    cdef double rot_mtrx[3][3]
    #cdef double [:,:] rot_mtrx = rot_mtrx_list # Can't use this line bc rot_mtrx needs to return a list
        
    cdef double M_anom, E_anom, T_anom, r_pl, r_star, mass_ratio
    
    cdef double r_vec_list[3]
    cdef double [:] r_vec = r_vec_list
    
    cdef double ang_pos[2]
    #cdef double [:] ang_pos = ang_pos_list
    
    #### POINT OF SPEED-UP: I could define C in an earlier function, where m*C = mass_ratio
    mass_ratio = ((m*M_jup)/(m_star*M_sun))
    
        
    M_anom = (2*pi/per)*elapsed_time
    
    # This is the eccentric anomaly at a given point in the epoch. It is different from the starting E_anomalies in E_anom_list in the make_arrays function.
    E_anom = kepler(M_anom, e)
    
    # T_anom replaces T_prog from the outdated code. It is the true anomaly after adding the randomly-sampled starting T_anom_0.
    T_anom = T_anom_0 + 2*atan( sqrt((1+e)/(1-e)) * tan(E_anom/2))

    rot_mtrx = rot_matrix(i, om, 0) # Omega = 0 arbitrarily 
    
    cdef double [:,:] rot_mtrx_memview = rot_mtrx # Convert list to memview here because matrix is all zeros if rot_matrix() returns a memview instead of a list. This is a costly line, but performin mat_mul below with meviews is ~5x as fast as with lists.
    

    ################### Angular Positions ######################
    
    r_pl = r(T_anom, a*au, e) # a is planet semi-major axis, so use it to find pl separation from bary with T_anom
    
    r_star = r_pl*mass_ratio
    

    # r_vec points from barycenter to the *star* (note the - sign) in the orbital plane, and has magnitude r_star. Like r_star, it has units of cm.
    r_vec[0] = -cos(T_anom)*r_star
    r_vec[1] = -sin(T_anom)*r_star
    r_vec[2] = 0
    
    # rot_vec points from barycenter to the star, but in coordinates where the xy-plane is the sky plane and the z-axis points toward Earth
    rot_vec = mat_mul(rot_mtrx_memview, r_vec)

    # ang_pos is the angular separation of the star from barycenter in milli-arcseconds.
    ang_pos[0] = rot_vec[0]*cm_2_mas
    ang_pos[1] = rot_vec[1]*cm_2_mas
    
    ################### Angular Velocities ########################
    
    #### POINT OF SPEED-UP: I don't need angular velocities for Hip, only ang_pos. So I could add an if statement here to pass if we're calculating for Hip, saving all the calculations below. This could give a ~25% speed-up.
    
    cdef double v_vec_pl[3]
    cdef double v_vec_star_list[3]
    cdef double [:] v_vec_star = v_vec_star_list
    
    cdef double rotated_v_vec_list[3]
    cdef double mu[2]
    
    v_vec_pl = v_vec(a, per, e, T_anom)
    # Stellar velocity is related to planet through their masses. 
    # Also they are in opposite directions, so add negative, but it shouldn't affect the final answer.
    v_vec_star[0] = -v_vec_pl[0] * mass_ratio
    v_vec_star[1] = -v_vec_pl[1] * mass_ratio
    v_vec_star[2] = -v_vec_pl[2] * mass_ratio
    
    
    rotated_v_vec = mat_mul(rot_mtrx_memview, v_vec_star)
    
    # mu is the proper motion of the star due to the planet's orbit in milli-arcseconds per year.
    mu[0] = rotated_v_vec[0]*cmd_2_masyr
    mu[1] = rotated_v_vec[1]*cmd_2_masyr

    return ang_pos, mu


def post_tot(double [:] rv_post_list, double [:] astro_post_list, int grid_num,
            long [:] a_inds, long [:] m_inds):
    
    cdef int size, i
    cdef double [:,:] tot_prob_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)
    cdef double prob
    
    size = rv_post_list.size
    
    for i in range(size):
    
        a_i = a_inds[i]
        m_i = m_inds[i]
    
        prob = rv_post_list[i]*astro_post_list[i]
    
        tot_prob_array[m_i, a_i] += prob
    
    return tot_prob_array

#@profile
cdef rot_matrix(double i, double om, double Om):
    """
    This is P3*P2*P1 from Murray & Dermott. It is not given explicitly in the text. They multiply it immediately by r*[cos(f), sin(f), 0]
    because this gives the projection of position onto the sky. However, we also need the projection of velocity, so we need the matrix
    pre-multiplication by the position vector.
    """
    cdef double sin_Om, sin_om, sin_i, cos_Om, cos_om, cos_i
    
    cdef double list_1[3]
    cdef double list_2[3]
    cdef double list_3[3]

    cdef double rot_mtrx[3][3] # This is a list
    #cdef double [:,:] rot_mtrx_mem = rot_mtrx # This line converts our list into a bytes-like object
    
    
    sin_Om = sin(Om)
    sin_om = sin(om)
    sin_i = sin(i)
    cos_Om = cos(Om)
    cos_om = cos(om)
    cos_i = cos(i)
    
    list_1 = [cos_Om*cos_om - sin_Om*cos_i*sin_om,
            -sin_om*cos_Om - sin_Om*cos_i*cos_om,
             sin_Om*sin_i]
    
    list_2 = [sin_Om*cos_om + cos_Om*cos_i*sin_om,
            -sin_om*sin_Om + cos_Om*cos_i*cos_om,
            -cos_Om*sin_i]
    
    list_3 = [sin_i*sin_om, sin_i*cos_om, cos_i]

    
    rot_mtrx[0][:] = list_1
    rot_mtrx[1][:] = list_2
    rot_mtrx[2][:] = list_3
    
    return rot_mtrx # Return a list, not a memview


cdef double r(double nu, double a, double e):
    """

    Equation of an ellipse (Murray & Dermott equation 2.20).
    Arguments:

        nu (radians): True anomaly
        a (distance): Semi-major axis of ellipse. Choice of a determines what output r represents. 
                        For example, if a is the semi-major axis of one planet's orbit, then r represents 
                        that planet's distance from barycenter as a function of nu. On the other hand,
                        if a is the SA of the test mass μ's orbit, then r is r1+r2 as a function of nu, 
                        where r1 (r2) is the distance of m1 (m2) from the system barycenter in the 
                        2-body (m1 & m2) frame.
        e (unitless): Eccentricity

    returns:
        r (same as a): Distance of particle from barycenter along its orbit
    """
    cdef double num, denom

    num = a*(1-e**2)
    denom = 1 + e*cos(nu)

    return num/denom
    
cdef mat_mul(double [:,:] mat, double [:] vec):
    """
    This is written specifically to matrix multiply rot_matrix (3x3) with 
    r_unit_vec (3x1) and v_vec_star (3x1) in astro_post_dense_loop.
    This function returns a list.
    """

    cdef int i, k
    
    cdef double result[3]     
                
    for i in range(3):
        result[i] = 0
        for k in range(3):
            result[i] += mat[i][k]*vec[k]

    return result

cdef v_vec(double a, double per, double e, double nu):
    """
    Uses Murray & Dermott equation 2.36. r_dot is not what we want because it doesn't capture the velocity perpendicular to the radial vector.
    Instead, v is the total velocity of the object. M&D doesn't actually give v vector explicitly, but I believe it's v_vec = [x_dot, y_dot, 0].

    Since periods created in units of days, v_vec has units of cm/day.
    v_vec is a list.
    """
    cdef double n, n_a, e_term, x_dot, y_dot
    
    cdef double v_vec[3]
    
    
    n = 2*pi/per
    n_a = n*a
    e_term = sqrt(1-e**2)

    x_dot = -n_a / e_term * sin(nu)
    y_dot = +n_a / e_term * (e + cos(nu))

    # To get the proper shape vector at the end, we need our 0 element to be an array with matching shape
    # zero_shape = np.shape(y_dot)
    # v_vec = np.array([x_dot, y_dot, np.zeros(zero_shape)])
    
    v_vec[0] = x_dot
    v_vec[1] = y_dot
    v_vec[2] = 0

    return v_vec

def r_dot(nu, a, P, e):
    """
    Murray & Dermott equation 2.31. This function gives the time rate of change
    of the distance between an orbiting body and the center of mass as a function of the body's true anomaly nu, period P, and
    eccentricity e.
    """
    num = 2*np.pi*a*e*np.sin(nu)
    denom = P*np.sqrt(1-e**2)

    return num/denom



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
    Given a 2D probability array, this function collapses the array along each axis to find the 68% confidence interval.

    value_spaces represents the parameter intervals covered by the array along each axis.
    It is expected in the form [(min_value1, max_value1), (min_value2, max_value2)], where 1 and 2 refer to the 0th and 1st axes.
    Note that the limits MUST be in this order: if the array has shape (x_num, y_num), then value_spaces must be [x_lims, y_lims].
    """
    bounds_list = []
    for i in range(2):

        array_1D = prob_array.sum(axis=i)
        grid_num = len(array_1D)


        # This gives only the 2-sigma, so that we get the 2-sigma limits at the end
        sig2 = contour_levels_1D(array_1D, [2])[0]

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
        
        
        inds_2sig = np.where(abs(interp_vals - sig2) < 1e-2*sig2)[0].take((0,-1))

        # value_bounds is a tuple of actual values, not indices
        value_bounds = index2value(inds_2sig, (0, interp_num-1), value_spaces[::-1][i])

        bounds_list.append(value_bounds)

    return bounds_list


def value2index(value, index_space, value_space):
    """
    The inverse of index2value: take a value on a 
    log scale and convert it to an index. index_space 
    and value_space are expected as tuples of the form 
    (min_value, max_value).
    """

    value = np.array(value)

    min_index, max_index = index_space[0],  index_space[1]
    min_value, max_value = value_space[0], value_space[1]

    index_range = max_index - min_index
    log_value_range = np.log10(max_value) - np.log10(min_value)

    index = (np.log10(value)-np.log10(min_value))*(index_range/log_value_range) + min_index

    return int(np.around(index))

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

    min_index, max_index = index_space[0],  index_space[1]
    min_value, max_value = value_space[0], value_space[1]

    index_range = max_index - min_index
    log_value_range = np.log10(max_value) - np.log10(min_value)

    # Convert from a linear space of indices to a linear space of log(values).
    log_value = (index-min_index)*(log_value_range/index_range) + np.log10(min_value)

    value = np.around(10**(log_value), 2) # Round to 2 decimal places

    return value

##################################################################################################


#@profile
#@cython.cdivision(True)
#@cython.boundscheck(False)
def gamma_direct(np.ndarray[double, ndim=1] a, np.ndarray[double, ndim=1] Mp, 
                 np.ndarray[double, ndim=1] per, np.ndarray[double, ndim=1] e, 
                 np.ndarray[double, ndim=1] i, np.ndarray[double, ndim=1] om,
                 np.ndarray[double, ndim=1] E):
    """
    Function to analytically calculate the first and second derivatives of the RV curve at a given point in the orbit.
    All arguments can be given as arrays (of compatible dimensions). M_anom and e in particular must be lists to work with the
    C-based kepler solver.
    Mp is expected in Jupiter masses.
    a is expected in au
    per is expected in days

    Returns:
    gamma_dot (m/s/d)
    gamma_ddot (m/s/d^2)
    """
    print('Using gamma_direct')
    cdef int size, j, k

    size = a.shape[0]

    cdef np.ndarray[double, ndim=1] gamma_dot = np.ndarray(shape=(size,), dtype=np.float64),\
                                    gamma_ddot = np.ndarray(shape=(size,), dtype=np.float64)
                                    
    cdef double cms2msday, cos_E, cos_nu, sin_nu_om, sin_E
    
    for k in range(size):
        Mp[k] = Mp[k]*M_jup
        a[k] = a[k]*au
    
    cms2msday = (1/100) * (24*3600)
    
    
    cdef double nu, E_dot, nu_dot, prefac, gd_t1, gd_t2, gd_t1_dot, gd_t2_dot, gdd_t1, gdd_t2
    
    for j in range(size):
        nu = 2*atan(((1+e[j])/(1-e[j]))**0.5*tan(E[j]/2))
        
        cos_E = cos(E[j])
        cos_nu = cos(nu)
        sin_nu_om = sin(nu+om[j])
        sin_E = sin(E[j])

        # Differentiate Kepler's equation in time to get E_dot
        # Note that E_dot has units of (1/per), where [per] is days. Therefore [gamma_ddot] = m/s/d^2
        E_dot = (2*pi/per[j])/(1-e[j]*cos_E)
        nu_dot = (1+tan(nu/2)**2)**-1 * ((1+e[j])/(1-e[j]))**0.5 * cos(E[j]/2)**-2 * E_dot

        # Convert prefac units from cm/s^2 to m/s/day
        # Negative just depends on choice of reference direction. I am being consistent with radvel rv_drive function.
        prefac = -(Mp[j]*G*sin(i[j]))/(a[j]**2*(1-e[j])) * cms2msday


        gd_t1 = (1+cos_nu)/(1+cos_E)
        gd_t2 = sin_nu_om/(1-e[j]*cos_E)


        gamma_dot[j] = prefac*gd_t1*gd_t2

        gd_t1_dot = ((1+cos_nu)*sin_E * E_dot - (1+cos_E)*sin(nu)*nu_dot) / (1+cos_E)**2
        gd_t2_dot = ((1-e[j]*cos_E)*cos(nu+om[j]) * nu_dot - sin_nu_om*e[j]*sin_E*E_dot) / (1-e[j]*cos_E)**2


        gdd_t1 = gd_t2 * gd_t1_dot
        gdd_t2 = gd_t1 * gd_t2_dot

        gamma_ddot[j] = prefac*(gdd_t1+gdd_t2)
    

    return gamma_dot, gamma_ddot

#@profile
@cython.cdivision(True)
@cython.boundscheck(False)
def gamma_direct_FAST(np.ndarray[double, ndim=1] a, np.ndarray[double, ndim=1] Mp, 
                      np.ndarray[double, ndim=1] per, np.ndarray[double, ndim=1] e, 
                      np.ndarray[double, ndim=1] i, np.ndarray[double, ndim=1] om, 
                      np.ndarray[double, ndim=1] E):
    """
    Function to analytically calculate the first and second derivatives of the RV 
    curve at a given point in the orbit. All arguments can be given as arrays 
    (of compatible dimensions). M_anom and e in particular must be lists to work with the
    C-based kepler solver.
    Mp is expected in Jupiter masses.
    a is expected in au
    per is expected in days

    Returns:
    gamma_dot (m/s/d)
    gamma_ddot (m/s/d^2)
    """
    print('Using gamma_direct_FAST')
    cdef int size

    size = a.shape[0]


    cdef np.ndarray[double, ndim=1] gamma_dot  = np.ndarray(shape=(size,), dtype=np.float64),\
                                    gamma_ddot = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #nu         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #nu_dot     = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #cos_E         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #cos_nu         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #sin_nu_om         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #sin_E         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #E_dot         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #prefac         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #gd_t1         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #gd_t2         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #gd_t1_dot         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #gd_t2_dot         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #gdd_t1         = np.ndarray(shape=(size,), dtype=np.float64),\
                                    #gdd_t2         = np.ndarray(shape=(size,), dtype=np.float64)

    Mp = Mp*M_jup
    a = a*au


    nu = 2*np.arctan(((1+e)/(1-e))**0.5*np.tan(E/2))

    cos_E = np.cos(E)
    cos_nu = np.cos(nu)
    sin_nu_om = np.sin(nu+om)
    sin_E = np.sin(E)

    # Differentiate Kepler's equation in time to get E_dot
    # Note that E_dot has units of (1/per), where [per] is days. Therefore [gamma_ddot] = m/s/d^2
    E_dot = (2*pi/per)/(1-e*cos_E)
    nu_dot = (1+np.tan(nu/2)**2)**-1 * ((1+e)/(1-e))**0.5 * np.cos(E/2)**-2 * E_dot

    # Convert prefac units from cm/s^2 to m/s/day
    # Negative just depends on choice of reference direction. I am being consistent with radvel rv_drive function.
    prefac = -(Mp*G*np.sin(i))/(a**2*(1-e)) * 864 # Save time on (24*3600)/100


    gd_t1 = (1+cos_nu)/(1+cos_E)
    gd_t2 = sin_nu_om/(1-e*cos_E)


    gamma_dot = prefac*gd_t1*gd_t2

    gd_t1_dot = ((1+cos_nu)*sin_E * E_dot - (1+cos_E)*np.sin(nu)*nu_dot) / (1+cos_E)**2
    gd_t2_dot = ((1-e*cos_E)*np.cos(nu+om) * nu_dot - sin_nu_om*e*sin_E*E_dot) / (1-e*cos_E)**2


    gdd_t1 = gd_t2 * gd_t1_dot
    gdd_t2 = gd_t1 * gd_t2_dot

    gamma_ddot = prefac*(gdd_t1+gdd_t2)


    return gamma_dot, gamma_ddot









