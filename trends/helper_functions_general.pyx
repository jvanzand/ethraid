# cython: language_level=3, boundscheck=False, cdivision=True, wraparound=False
## cython: binding=True

import numpy as np
cimport numpy as np
import scipy as sp
import scipy.stats as spst

import cython
from libc.math cimport sin, cos, tan, atan, sqrt, log

cdef float pi, two_pi, math_e, G, M_sun, M_jup, au, pc_in_cm

pi = 3.141592653589793
two_pi = 6.283185307179586
math_e = 2.718281828459045
G = 6.674299999999999e-08
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
au = 14959787070000.0
pc_in_cm = 3.086e18

np.random.seed(0)
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

    #np.random.seed(0)
    tp = 0
    a_min = a_lim[0]
    a_max = a_lim[1]
    m_min = m_lim[0]
    m_max = m_lim[1]


    # These semimajor axes are distances between the planet and the barycenter of the system. The star is on its own orbit, which we will get later.
    a_list = spst.loguniform.rvs(a_min, a_max, size=num_points)
    m_list = spst.loguniform.rvs(m_min, m_max, size=num_points)

    # Match up a_list and m_list and get the period for each pair (in days).
    per_list = P_list(a_list, m_list, m_star) # Use this line when we are actually sampling a_tot, not a_planet
    #per_list = P(a_list * (m_star+m_list*(M_jup/M_sun))/m_star, m_star+m_list*(M_jup/M_sun) )
    
    # Eccentricities drawn from a beta distribution. I am using (a,b) = (0.867, 3.03) according to Winn & Fabrycky (2014).
    e_list = spst.beta(0.867, 3.03).rvs(num_points)
    e_list = np.where(e_list > 0.99, 0.99, e_list) # Replace e > 0.99 with 0.99

    cosi_list = np.random.uniform(0, 1, num_points)
    i_list = np.arccos(cosi_list)

    # Mean anomaly, uniformly distributed. This represents M at the beginning of the Hipparcos epoch for BOTH RVs and astrometry. Use this to solve for True anomaly.
    M_anom_0_list = np.random.uniform(0, two_pi, num_points)

    # Evolving M_anom forward to the epoch of RV calculations.
    # Commenting out because I evolve forward individually in the RV helper module
    #M_anom_evolved = M_anom_0 + two_pi*((rv_epoch - hip_times[0])/per_list)
    #E_anom_rv = ck.kepler_array(M_anom_evolved, e_list) # Used in post_rv

    # Arguments of peri, uniformly distributed
    om_list = np.random.uniform(0, two_pi, num_points)

    # Breaking up the (a, M) parameter space into grid_num x grid_num
    a_bins = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
    m_bins = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)

    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)

    return a_list, m_list, per_list, e_list, i_list, om_list, M_anom_0_list, a_inds, m_inds


def post_tot(double [:] rv_post_list, double [:] astro_post_list, int grid_num,
            long [:] a_inds, long [:] m_inds):
            
    """
    Start with 2 1D lists and multiply them element-wise. 
    THEN form the result into a 2D array.
    """

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


def P_list(double [:] a_list, double [:] m_list, double m_star):
    """
    Uses Kepler's third law to find the periods (in days) of a 
    list of planet masses (in M_jup) given their semimajor axes (au)
    and the mass of a star (M_sun, note that m_star is a float, not a list).
    
    a (au): list of semi-major axes
    Mp (M_Jup): list of companion masses
    Ms (M_sun): a single stellar mass
    """

    cdef int length, j
    cdef double a, m
    length = a_list.size
    
    cdef np.ndarray[double, ndim=1] per_list = np.ndarray(shape=(length,), dtype=np.float64)
    
    for j in range(length):
        a = a_list[j]
        m = m_list[j]
        per_list[j] = P(a, m, m_star)
        
    return per_list

cpdef P(double a, double m, double m_star):
    """
    Uses Kepler's third law to find the period of a planet (in days) given its
    semimajor axis, the planet mass, and the stellar mass.

    a (au): semi-major axis
    Mp (M_Jup): companion mass
    Ms (M_sun): stellar mass
    """
    
    cdef double m_g, m_star_g, sec_2_days, P_days

    m_g = m*M_jup
    m_star_g = m_star*M_sun

    sec_2_days = 1./(24*3600) # Note the 1.; with 1, the result would be 0

    P_days = sqrt((two_pi)**2*(a*au)**3/(G*(m_g + m_star_g))) * sec_2_days

    return P_days
    

def contour_levels(prob_array, sig_list, t_num = 1e3):
    """
    Contour drawing method taken from 
    https://stackoverflow.com/questions/37890550/python-plotting-percentile-contour-lines-of-a-probability-distribution
    This function takes a 2-D array of probabilities and returns a 1-D array 
    of the probability values corresponding to 1-sigma and 2-sigma contours. 
    In this case, the 1-sigma contour encloses 68% of the total probability. 
    The array is expected to be normalized. sig_list is a list containing 
    any combination of the integers 1, 2, or 3 to indicate desired contours. 
    For example, [1,3] will return the 1 and 3 sigma contours.
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
    Same as contour_levels, but adapted for 1D arrays. 
    Hopefully I can condense these into 1 in the future.
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

    index_range = max_index - min_index
    log_value_range = np.log10(max_value) - np.log10(min_value)

    value_arr = np.array(value)

    index = (np.log10(value_arr)-np.log10(min_value))\
                                    *(index_range/log_value_range) + min_index

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

    min_index, max_index = index_space[0],  index_space[1]
    min_value, max_value = value_space[0], value_space[1]

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


    return a_cm / au
    
    