from tqdm import tqdm
import numpy as np
from astroquery.vizier import Vizier

import cython
cimport numpy as np
from libc.math cimport sin, cos, tan, atan, sqrt, log

from ethraid import hip_times, hip_epoch, gaia_times, gaia_epoch
from ethraid.compiled._kepler import kepler_single

cdef double two_pi, math_e, G, M_sun, M_jup, au, pc_in_cm, baseline_yrs

pi = 3.141592653589793
two_pi = 6.283185307179586
math_e = 2.718281828459045

# G in AU, M_Jup, day units.
G = 2.824760877012879e-07 # (c.G.cgs*(1/c.au.cgs)**3 * (c.M_jup.cgs) * (24*3600)**2).value

baseline_yrs = gaia_epoch - hip_epoch

def astro_list(double [:] a_list, double [:] m_list, double [:] e_list, 
               double [:] i_list, double [:] om_list, double [:] M_anom_0_list, 
               double [:] per_list, 
               double m_star, double d_star, double delta_mu, double delta_mu_err):
    """
    Calculates the log-likelihood of the astrometry data conditioned 
    on each of a list of orbital models, resulting in a list of 
    log-likelihoods. All input lists must have the same length, 
    which is also the length of the output, log_lik_list.
    
    Arguments:
        a_list (list of floats, AU): Semi-major axes
        m_list (list of floats, M_jup): Companion masses
        e_list (list of floats): Orbital eccentricities
        i_list (list of floats, radians): Orbital inclinations
        om_list (list of floats, radians): Arguments of periastron
        M_anom_0_list (list of floats, radians): Mean anomalies at the 
                                                 beginning of the Hipparcos 
                                                 mission (1989.85)
        per_list (list of floats, days): Orbital periods
        m_star (float, M_jup): Host star mass
        d_star (float, AU): Distance of system from Earth
        delta_mu (float, mas/yr): Magnitude of difference 
                                  between measured Gaia pm and 
                                  positional average pm between 
                                  Gaia and Hipparcos.
        delta_mu_err (float, mas/yr): Error on delta_mu
    
    Returns:
        log_lik_list (list of floats): List of log_likelihoods corresponding
                                       to the input models. 
    """
      
    cdef int num_points, j
    cdef double a, m, e, i, om, M_anom_0, per, log_lik
    num_points = a_list.shape[0]     
        
    cdef np.ndarray[double, ndim=1] log_lik_list = np.ndarray(shape=(num_points,), dtype=np.float64)
    
    if delta_mu < 0:
       raise Exception("delta_mu cannot be negative")
    
    print('Running astrometry models')
    for j in tqdm(range(num_points)):
       
        a = a_list[j]
        m = m_list[j]
        e = e_list[j]
        i = i_list[j]
        om = om_list[j] # Because om_list contains the *planet* \omegas, this om should technically be om_list[j]+pi. However, the dmu() function below is insensitive to this change, so we use om_list[j] itself for simplicity.
        M_anom_0 = M_anom_0_list[j]
        per = per_list[j]

        log_lik = log_lik_dmu(a, m, e, i, om, M_anom_0,
                              per, m_star, d_star,
                              delta_mu, delta_mu_err)

        log_lik_list[j] = log_lik
    
    return log_lik_list


def log_lik_dmu(double a, double m, double e, double i, double om, double M_anom_0, 
                double per, double m_star, double d_star,
                double dmu_data, double dmu_data_err):
    """
    Computes the log-likelihood of the astrometry data conditioned on
    a given orbital model (set of a, Mp, e, i, om, and M_anom_0).
    
    Arguments:
        a (float, AU): Semi-major axis
        m (float, M_jup): Companion mass
        e (float): Orbital eccentricity
        i (float, radians): Orbital inclination
        om (float, radians): Argument of periastron
        M_anom_0 (float, radians): Mean anomaly at the beginning
                                   of the Hipparcos mission (1989.85)
        per (float, days): Orbital period
        m_star (float, M_jup): Host star mass
        d_star (float, AU): Distance of system from Earth
        dmu_data (float, mas/yr): Magnitude of difference 
                                  between measured Gaia pm and 
                                  positional average pm between 
                                  Gaia and Hipparcos.
        dmu_data_err (float, mas/yr): Error on dmu_data
    
    Returns:
        log_likelihood (float): Natural logarithm of the
                                model likelihood.
    """
    cdef double dmu_model, log_likelihood
    
    dmu_model = dmu(a, m, e, i, om, M_anom_0, per, m_star, d_star)
    
    # Log of the exponential. Omit prefactor because it has no effect on shape of likelihood surface.
    log_likelihood = -(dmu_data - dmu_model)**2/(2*dmu_data_err**2)
    
    return log_likelihood


def dmu(double a, double m, double e, double i, double om, double M_anom_0, 
        double per, double m_star, double d_star):
    """
    Computes the change in proper motion, delta_mu, for 
    a set of model parameters
    
    Arguments:
        a (float, AU): Semi-major axis
        m (float, M_jup): Companion mass
        e (float): Orbital eccentricity
        i (float, radians): Orbital inclination
        om (float, radians): Argument of periastron
        M_anom_0 (float, radians): Mean anomaly at the beginning
                                   of the Hipparcos mission (1989.85)
        per (float, days): Orbital period
        m_star (float, M_jup): Host star mass
        d_star (float, AU): Distance of system from Earth
    
    Returns:
        dmu_model (float, mas/yr): Modeled magnitude of difference 
                                   between Gaia pm and positional 
                                   average pm between Gaia and 
                                   Hipparcos. 
        
    """
    cdef double mass_ratio, au_2_mas, aud_2_masyr
    cdef double mean_motion, a_star, M1, M2, E1, E2
    cdef double x_pos_avg, y_pos_avg, x_vel_avg, y_vel_avg
    cdef double cos_om, sin_om, cos_i
    
    # vec_list is a temporary variable to define vec
    # vec holds various values throughout dmu(). After each value
    # has served its purpose, it is overwritten so that only one
    # vector needs to be allocated, saving time.
    cdef double vec_list[3]
    cdef double [:] vec = vec_list
    
    cdef double time_endpoints[2][2]
    cdef double ang_pos_avg[2][2]
    cdef double mu_avg[2][2]
    cdef double mu_gaia[2]
    cdef double mu_hg[2]
    
    d_star = d_star/206264.80624548031 # Divide by (c.pc.cgs/c.au.cgs).value to get units of pc
    au_2_mas = 1e3/d_star # Conversion factor btwn au and milli-arcseconds, with d_star converted to pc
    aud_2_masyr = au_2_mas * 365.25 # au/day to milli-arcseconds/year conversion factor

    time_endpoints = [[hip_times[0], gaia_times[0]], 
                      [hip_times[1], gaia_times[1]]]

    mass_ratio = m/(m_star + m)
    a_star = a*mass_ratio
    mean_motion = two_pi/per
    
    #rot_matrix(i, om, rot_mtrx)
    cos_om = cos(om)
    sin_om = sin(om)
    cos_i = cos(i)
    
    for l in range(2): # Hipparcos or Gaia
        start_time = time_endpoints[0][l] - time_endpoints[0][0] # The "start time" of Hip or Gaia relative to the start of Hip. For Hip, start_time is 0 by definition. For Gaia, it is the time between Hip_start and Gaia_start
        end_time = time_endpoints[1][l] - time_endpoints[0][0] # End time relative to the start of Hip.

        ## Mean anomaly is the elapsed time times the mean motion, plus a randomly-sampled starting mean anomaly
        M1 = mean_motion*start_time + M_anom_0
        M2 = mean_motion*end_time + M_anom_0

        E1 = kepler_single(M1, e)
        E2 = kepler_single(M2, e)
        
        # Get avg. separation of the star in au, but pointing in the wrong direction. E1 and E2 correspond to companion direction, so we need a factor of -1 to flip direction (and multiplying the tuple by -1 in this line gives an empty tuple).
        x_pos_avg, y_pos_avg = pos_avg(a_star, mean_motion, e, 
                                       E1, E2, start_time, end_time)
        
        # Two transformations are performed in this step. First, x_pos_abg and y_pos_avg are rotated into the sky plane using the components of the rotation matrix from Murray and Dermott Eqs. 2.119-2.121      
        # Second, the leading negative signs point vec toward the *star*, not the companion.
        vec[0] = -(x_pos_avg*cos_om - y_pos_avg*sin_om)
        vec[1] = -(x_pos_avg*cos_i*sin_om + y_pos_avg*cos_i*cos_om)
        vec[2] = 0

        # Average angular position of the star relative to barycenter in milli-arcseconds.
        # ang_pos is in the plane of the sky, so no need to use vec's third component, which points out of the sky
        ang_pos_avg[l][0] = vec[0]*au_2_mas
        ang_pos_avg[l][1] = vec[1]*au_2_mas
        
        ################### Angular Velocities ########################
        # Only need Gaia velocities, so skip Hip velocities
        if l == 0:
            continue

        # Get average velocity of the star (au/day)
        x_vel_avg, y_vel_avg = vel_avg(a_star, mean_motion, e, 
                                       E1, E2, start_time, end_time)
        
        # Like position, rotate the velocity components and flip them to point to the star.
        vec[0] = -(x_vel_avg*cos_om - y_vel_avg*sin_om)
        vec[1] = -(x_vel_avg*cos_i*sin_om + y_vel_avg*cos_i*cos_om)
        vec[2] = 0


        # mu_avg is a 2x2 array. The top row stays empty because we skip Hip. (The l==0 case is never executed
        # because we skip it in the if statement above.)
        # The bottom row is average proper motion over the Gaia mission.
        # The average proper motion of the star due to the companion's orbit is in milli-arcseconds per year.
        # Since period is in days and a_star_units is in au, velocities are in au/day. Convert to mas/yr.
        mu_avg[l][0] = vec[0]*aud_2_masyr
        mu_avg[l][1] = vec[1]*aud_2_masyr
        ###############################################################
        ###############################################################

    mu_gaia = mu_avg[1]
    # To get the positional avg., subtract the epoch positions and divide by the time between in years.
    # First index tells Hip ([0]) or Gaia ([1]), second index tells x ([0]) or y ([1])
    # Units of mas/yr
    mu_hg[0] = (ang_pos_avg[1][0] - ang_pos_avg[0][0])/baseline_yrs # x-comp. = gaia_x - hip_x
    mu_hg[1] = (ang_pos_avg[1][1] - ang_pos_avg[0][1])/baseline_yrs # y-comp. = gaia_y - hip_y
    
    dmu_model = sqrt((mu_hg[0] - mu_gaia[0])**2 + (mu_hg[1] - mu_gaia[1])**2) # Pythagorean distance
    
    return dmu_model

cdef pos_avg(double a, double n, double e, double E1, double E2, 
             double t1, double t2):
    """
    Calculate the average x/y positions of an object on an elliptical orbit, 
    where (0,0) is the focus and E = 0 points to +x. Equations for average
    position come from integrating Murray & Dermott Equations 2.41 and
    dividing by the time interval.

    Arguments:
        a (float, au): semi-major axis of the ellipse. Note that for an orbiting body,
                       a must be the semi-major axis of that body's orbital path,
                       not a_total of the full two-body orbit.
        n (float, radians per day): 2pi/per, aka mean motion
        e (float): orbital eccentricity
        E1, E2 (float, radians): eccentric anomaly at t1 and t2
        t1, t2 (float, days): beginning and ending time to calculate average

    Returns:
        x_avg, y_avg (float, au): average x and y positions
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
    where (0,0) is the focus and E = 0 points to +x.
    
    Arguments:
        a (float, au): semi-major axis of the ellipse. Note that for an orbiting body,
                       a must be the semi-major axis of that body's orbital path,
                       not a_total of the full two-body orbit.
        n (float, radians per day): 2pi/per, aka mean motion
        e (float): orbital eccentricity
        E1, E2 (float, radians): eccentric anomaly at t1 and t2
        t1, t2 (float, days): beginning and ending time to calculate average

    Returns: 
        x_avg, y_avg (float, au/day): average x and y velocities
    """

    cdef double x_term_1, x_term_2, x_integral, x_avg,\
                y_term_1, y_term_2, y_integral, y_avg


    x_term_1 = a*cos(E1)
    x_term_2 = a*cos(E2)
    x_integral = x_term_2 - x_term_1

    x_avg = 1/(t2-t1) * x_integral

    y_term_1 = a*sqrt(1-e**2)*sin(E1)
    y_term_2 = a*sqrt(1-e**2)*sin(E2)
    y_integral = y_term_2 - y_term_1

    y_avg = 1/(t2-t1) * y_integral


    return x_avg, y_avg

cdef void rot_matrix(double i, double om, double [:,:] rot_mtrx):
    """
    This is P2*P1 from Murray & Dermott Eqs. 2.119. We omit P3 (2.120) because
    the longitude of the ascending node (Omega) can be set arbitrarily to 0 for
    our purposes, saving some calculations.

    This function doesn't return anything. Instead, declare a matrix in your 
    function and this will update it, more time by not allocating memory to
    and returning a matrix.

    Arguments:
        i (float, radians): orbital inclination (0 = face-on)
        om (float, radians): argument of periastron
        rot_mtrx (3x3 array of zeros): Initial array, modified in place

    Returns:
        None (but populates rot_mtrx with new values)
    """
    cdef double sin_om, sin_i, cos_om, cos_i

    sin_om = sin(om)
    sin_i  = sin(i)
    cos_om = cos(om)
    cos_i  = cos(i)

    rot_mtrx[0][0] = cos_om
    rot_mtrx[0][1] = -sin_om
    rot_mtrx[0][2] = 0

    rot_mtrx[1][0] = cos_i*sin_om
    rot_mtrx[1][1] = cos_i*cos_om
    rot_mtrx[1][2] = -sin_i

    rot_mtrx[2][0] = sin_i*sin_om
    rot_mtrx[2][1] = sin_i*cos_om
    rot_mtrx[2][2] = cos_i

    # Do not return rot_mtrx
            

# THIS version of the function will let me input the same vector as in_vec and out_vec.
# I want to do this because it will allow me to allocate fewer arrays in the dmu function,
# offering speedups.
cpdef void mat_mul(double [:,:] mat, double [:] in_vec, double [:] out_vec):
    """
    This function is written specifically to matrix multiply rot_mtrx (3x3)
    with vec (3x1) in the dmu() function. Saves time by receiving its 
    "output" (out_vec) as an argument and modifying it in place without
    returning anything. I made this function even more ad-hoc by removing
    the for-loop over vector elements. This lets the user pass the same
    vector as both the input and the output, which saves time by only
    allocating one array.
    
    Arguments:
        mat_mul (3x3 array of floats): Rotation matrix from the rot_matrix
                                       function
        in_vec (3x1 array of floats): Position vector of star [x, y, z]
        
    Returns:
        None (but populates out_vec)
    """

    cdef double m00, m01, m02,\
                m10, m11, m12,\
                m20, m21, m22,\
                iv0, iv1, iv2,\
                ov0, ov1, ov2
    
    m00 = mat[0][0]
    m01 = mat[0][1]
    m02 = mat[0][2]
    m10 = mat[1][0]
    m11 = mat[1][1]
    m12 = mat[1][2]
    m20 = mat[2][0]
    m21 = mat[2][1]
    m22 = mat[2][2]
    
    iv0 = in_vec[0]
    iv1 = in_vec[1]
    iv2 = in_vec[2]
    
    ov0 = m00*iv0 + m01*iv1 + m02*iv2
    ov1 = m10*iv0 + m11*iv1 + m12*iv2
    ov2 = m20*iv0 + m21*iv1 + m22*iv2
    
    out_vec[0] = ov0
    out_vec[1] = ov1
    out_vec[2] = ov2
            
    # Do not return out_vec

def HGCA_retrieval(hip_id=None, gaia_id=None):
    """
    Access the HGCA (EDR3 version) on Vizier and calculate
    Δμ and Δμ_error for a given target. This function allows
    a user to input astrometry constraints without having
    to download the HGCA or propagate errors themselves.
    
    NOTE: Vizier's values are rounded to 3 decimal places. 
          The Vizier quantities are still consistent with
          those from the HGCA, but not perfectly equal.
          
    Arguments:
        hip_id (str): Hipparcos identifier
        gaia_id (str): Gaia identifier
    
    Returns:
        dmu (float, mas/yr): Absolute value of the difference between the 
                             average Hipparcos-Gaia proper motion over the 
                             ~24 years between the two missions and the 
                             Gaia proper motion.
        dmu_err (float, mas/yr): Error on dmu
    """
    if hip_id is None and gaia_id is None:
        print("helper_functions_astro.HGCA_retrieval: No HGCA identifier given.")
        dmu = 0
        dmu_err = 1e8
        return dmu, dmu_err

    if hip_id is not None and gaia_id is not None:
        filter_dict = {'HIP':hip_id, 'Gaia':gaia_id}
    elif hip_id is not None:
        filter_dict = {'HIP':hip_id}
    elif gaia_id is not None:
        filter_dict = {'Gaia':gaia_id}

    # pmRA and pmDE are the Gaia proper motions 
    include_cols = ['HIP', 'Gaia',
                    'pmRA',  'pmDE',   'pmRAhg',   'pmDEhg', 
                    'e_pmRA','e_pmDE', 'e_pmRAhg', 'e_pmDEhg',
                    'EpochRAgaia', 'EpochDEgaia', 'EpochRAhip', 'EpochDEhip']

    # Query Vizier for HGCA entry. 'J/ApJS/254/42' is the name of the
    # HGCA EDR3.
    v = Vizier(columns=include_cols,
               column_filters=filter_dict, catalog='J/ApJS/254/42')
    
    table_set = v.get_catalogs(v.catalog)
    
    if len(table_set)==0:
        raise Exception("helper_functions_astro.HGCA_retrieval:\n"
                        "           No matching targets found. Target may not be in the HGCA.")

    table = v.get_catalogs(v.catalog)[0]
    
    if len(table)>1:
        raise Exception("helper_functions_astro.HGCA_retrieval:\n"
                        "           Multiple matching targets found in HGCA.\n"
                        "           Try a different identifier or enter Δμ manually.")
    
    pmra_gaia = table['pmRA']
    pmra_hg = table['pmRAhg']
    pmdec_gaia = table['pmDE']
    pmdec_hg = table['pmDEhg']

    pmra_gaia_error = table['e_pmRA']
    pmra_hg_error = table['e_pmRAhg']
    pmdec_gaia_error = table['e_pmDE']
    pmdec_hg_error = table['e_pmDEhg']


    dmu = float(calc_dmu(pmra_gaia, pmra_hg, pmdec_gaia, pmdec_hg))
    dmu_err = float(calc_dmu_error(pmra_gaia_error, pmra_hg_error,\
                                   pmdec_gaia_error, pmdec_hg_error,\
                                   pmra_gaia, pmra_hg, pmdec_gaia, pmdec_hg))
    
    return dmu, dmu_err
    

def calc_dmu(pmra_gaia, pmra_hg, pmdec_gaia, pmdec_hg):
    """
    Simple pythagorean calculation of |proper motion change|
    
    Arguments:
        pmra_gaia (float, mas/yr): Proper motion in RA during the Gaia mission
        pmra_hg (float, mas/yr): Average proper motion in RA between the 
                                 Hipparcos and Gaia missions
        pmdec_gaia (float, mas/yr): Proper motion in DEC during the Gaia mission
        pmdec_hg (float, mas/yr): Average proper motion in DEC between the 
                                  Hipparcos and Gaia missions
    
    Returns:
        dmu_mag (float, mas/yr): Absolute value of the difference between the 
                                 average Hipparcos-Gaia proper motion over the 
                                 ~24 years between the two missions and the 
                                 Gaia proper motion.
    """
    
    dmu_mag = np.sqrt((pmra_gaia-pmra_hg)**2+(pmdec_gaia-pmdec_hg)**2)

    return dmu_mag

def calc_dmu_error(pmra_gaia_err, pmra_hg_err, pmdec_gaia_err, pmdec_hg_err, pmra_gaia, pmra_hg, pmdec_gaia, pmdec_hg):
    """
    Error on proper motion magnitude using error propagation formula.
    This simplified formula assumes independent variables.
    
    Arguments:
        pmra_gaia_err (float, mas/yr): Error on pmra_gaia
        .
        .
        .
    
    Returns:
        dmu_err (float, mas/yr): Error on the change in proper motion
    """

    numerator = np.sqrt((pmra_gaia-pmra_hg)**2   * (pmra_gaia_err**2+pmra_hg_err**2)\
                      + (pmdec_gaia-pmdec_hg)**2 * (pmdec_gaia_err**2+pmdec_hg_err**2))

    denominator = calc_dmu(pmra_gaia, pmra_hg, pmdec_gaia, pmdec_hg)
    
    dmu_err = numerator/denominator

    return dmu_err





