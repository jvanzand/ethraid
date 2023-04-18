import sys
import numpy as np
import scipy as sp
import scipy.stats as spst

import radvel as rv
cimport numpy as np
cimport cython
from libc.math cimport sin, cos, tan, atan, sqrt, log

cdef double pi, two_pi, math_e, G, M_sun, M_jup, au, pc_in_cm

pi = 3.141592653589793
two_pi = 6.283185307179586
math_e = 2.718281828459045
old_G = 6.674299999999999e-08

G = 2.824760877012879e-07 # c.G.cgs.value*(1/c.au.cgs.value)**3 * (c.M_jup.cgs.value) * (24*3600)**2

M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
au = 14959787070000.0
pc_in_cm = 3.086e18

def make_arrays(double m_star, tuple a_lim, tuple m_lim, int grid_num, int num_points):
    """
    Create the parameter arrays which will be used for the RV and astrometry posteriors.
    
    Arguments:
        m_star (float, M_sun): Mass of host star
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max)
        m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
        grid_num (int): Dimensions of square (a,m) grid
        num_points (int): Number of random orbital models to simulate
    
    Returns:
        a_list, m_list, per_list, e_list, 
        i_list, om_list, M_anom_0_list (numpy arrays, len = num_points):
                                        Lists of randomly sampled semi-major axis, mass,
                                        eccentricity, inclination, argument of
                                        periastron, and initial mean anomaly. Do not
                                        sample directly in period.
        a_inds, m_inds (numpy arrays of ints, len = num_points): Grid position where each 
                                        (a, m, per, e, i, om, M_anom_0) model 
                                        will be placed, based on the model's 
                                        a and m values
    """

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

    np.random.seed(10)
    tp = 0
    a_min = a_lim[0]
    a_max = a_lim[1]
    m_min = m_lim[0]
    m_max = m_lim[1]


    # These are the "full" semi-major axes of the orbit, ie the sma of the ellipse traced by the 1-body solution to the 2-body problem. a = a_planet+a_star
    a_list = spst.loguniform.rvs(a_min, a_max, size=num_points)
    m_list = spst.loguniform.rvs(m_min, m_max, size=num_points)

    # Match up a_list and m_list and get the period for each pair (in days).
    # Calculate this now to avoid having to do it twice for RVs and astrometry.
    per_list = P_list(a_list, m_list, m_star) # Use this line when we are sampling a_tot, not a_planet
    
    # Eccentricities drawn from a beta distribution.
    #e_list = np.zeros(num_points)
    #e_list = spst.beta(0.867, 3.03).rvs(num_points) # Beta distribution for planets+BDs together from Bowler+2020
    #e_list = np.where(e_list > 0.99, 0.99, e_list) # Replace e > 0.99 with 0.99
    
    e_list = ecc_dist(per_list, num_points)

    cosi_list = np.random.uniform(0, 1, num_points)
    i_list = np.arccos(cosi_list)

    # Mean anomaly, uniformly distributed. This represents M at the beginning of the Hipparcos epoch for BOTH RVs and astrometry. Use this to solve for True anomaly.
    M_anom_0_list = np.random.uniform(0, two_pi, num_points)

    # Arguments of peri, uniformly distributed
    om_list = np.random.uniform(0, two_pi, num_points)

    # Breaking up the (a, M) parameter space into grid_num x grid_num
    a_bins = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
    m_bins = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)

    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)

    return a_list, m_list, per_list, e_list, i_list,\
           om_list, M_anom_0_list, a_inds, m_inds

def post_tot(double [:] rv_post_list, double [:] astro_post_list, 
            double [:] imag_post_list, int grid_num,
            long [:] a_inds, long [:] m_inds):

   """
   Start with 2 1D lists and multiply them element-wise, THEN form 
   the result into a 2D array. This function is for the total posterior;
   the individual RV, astrometry, and imaging posteriors are handled by the 
   post_single() function below.

   Arguments:
       rv_post_list (np array of floats, len=num_points): List of model likelihoods
                    given the RV data
       rv_post_list (np array of floats, len=num_points): List of model likelihoods
                   given the astrometry data
       grid_num (int): Dimension of square (a,m) grid
       a_inds, m_inds (numpy arrays of ints, len = num_points): Grid position where each 
                                       (a, m, per, e, i, om, M_anom_0) model will 
                                       be placed, based on the model's 
                                       a and m values
                            
   Returns:
       tot_prob_array (numpy array, dim = grid_num x grid_dum): 2-D array of binned
                      (aka marginalized) posterior probabilities. Note, the binning
                      process itself is what applies my priors, converting the
                      individual likelihoods into posterior probabilities.
   """

   cdef int num_points, i, a_i, m_i
   cdef double prob
   
   # Not cdef because then I can't change it from memview to np array
   tot_prob_array = np.zeros((grid_num, grid_num))

   num_points = rv_post_list.size

   for i in range(num_points):

       a_i = a_inds[i]
       m_i = m_inds[i]

       prob = rv_post_list[i]*astro_post_list[i]*imag_post_list[i]

       tot_prob_array[m_i, a_i] += prob


   tot_prob_array = np.array(tot_prob_array)

   return tot_prob_array/tot_prob_array.sum()

def post_tot_simplified(double [:] rv_post_list, double [:] astro_post_list, 
                        double [:,:] post_imag, int grid_num,
                        long [:] a_inds, long [:] m_inds):
            
    """
    Start with 2 1D lists and multiply them element-wise, THEN form 
    the result into a 2D array. This function is for the total posterior;
    the individual RV, astrometry, and imaging posteriors are handled by the 
    post_single() function below.
    
    Arguments:
        rv_post_list (np array of floats, len=num_points): List of model likelihoods
                     given the RV data
        rv_post_list (np array of floats, len=num_points): List of model likelihoods
                    given the astrometry data
        grid_num (int): Dimension of square (a,m) grid
        a_inds, m_inds (numpy arrays of ints, len = num_points): Grid position where each 
                                        (a, m, per, e, i, om, M_anom_0) model will 
                                        be placed, based on the model's 
                                        a and m values
                                        
    Returns:
        tot_prob_array (numpy array, dim = grid_num x grid_dum): 2-D array of binned
                       (aka marginalized) posterior probabilities. Note, the binning
                       process itself is what applies my priors, converting the
                       individual likelihoods into posterior probabilities.
    """

    cdef int size, i, a_i, m_i
    cdef double [:,:] rv_ast_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)
    cdef double prob

    size = rv_post_list.size

    for i in range(size):

        a_i = a_inds[i]
        m_i = m_inds[i]

        prob = rv_post_list[i]*astro_post_list[i]

        rv_ast_array[m_i, a_i] += prob
    
    # Not cdef because then I can't change it from memview to np array
    tot_prob_array = np.zeros((grid_num, grid_num))
    
    # post_imag is not a list like RVs and astrometry above. It is input to this function as a 2D array. This is because it's a lot easier to calculate. It would take way longer if we calculated a length-1e6 (or 1e8) list instead of a 100x100 array.
    for i in range(grid_num):
        for j in range(grid_num):
            tot_prob_array[i,j] = rv_ast_array[i,j]*post_imag[i,j]
    
    
    tot_prob_array = np.array(tot_prob_array)

    return tot_prob_array/tot_prob_array.sum()


def post_single(double [:] prob_list, long [:] a_inds, long [:] m_inds, int grid_num):
    """
    Form a list of probabilities into a 2D array.
    
    Arguments:
        prob_list (list/array): List of probabilities to be be reshaped
        a_inds, m_inds (numpy arrays of ints): Coordinates where each probability
                                               will be placed. Probabilities with
                                               matching coordinates are summed 
                                               (marginalized).
        grid_num (int): Dimension of square array into which prob_list will be 
                        formed
    
    Returns:
        prob_array (np array, dim = grid_num x grid_num): Array of binned 
                                                         probabilities
    """

    cdef int i, size, a_i, m_i
    cdef double [:,:] prob_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)

    size = prob_list.shape[0]

    for i in range(size):

        a_i = a_inds[i]
        m_i = m_inds[i]

        prob_array[m_i, a_i] += prob_list[i]
        
    np_prob_array = np.array(prob_array)/np.array(prob_array).sum()
    
    return np_prob_array


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
        mp = m_list[j]
        per_list[j] = P(a, mp, m_star)
        
    return per_list

cpdef P(double a, double m_planet, double m_star):
    """
    Uses Kepler's third law to find the period of a planet (in days) given its
    semimajor axis, planet mass, and stellar mass.
    
    Arguments:
        a (float, au): Semi-major axis
        m_planet (float, M_jup): Planet mass
        m_star (float, M_jup): Stellar mass
    
    Returns:
        per (float, days): Companion orbital period
    """
    
    cdef double per
    
    per = two_pi * a**(1.5) * (G*(m_planet+m_star))**(-0.5)
    
    return per


def ecc_dist(double [:] per_list, int num_points):
    """
    Sample a random eccentricity whose distribution is based on a and m.
    Kipping(2013) advocates two distributions for P below and above 382.3 days.
    
    It might make more sense to use a third (Bowler) distribution 
    for much longer periods (bc Kipping only used 400 exoplanet 
    eccentricities back in 2013, and no BDs).
    
    Arguments:
        per_list (list of floats, days): List of companion periods
        num_points (int): Number of random eccentricities to generate
    
    Returns:
        e_list (list of floats): List of randomly sampled eccentricities
    """
    cdef int i
    cdef double per, e
    
    cdef np.ndarray[double, ndim=1] e_list = np.ndarray(shape=(num_points), dtype=np.float64),\
                                    kipping_short = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    kipping_long = np.ndarray(shape=(int(num_points)), dtype=np.float64)
    
          
    # Note that I make lists that are each num_points long, so I only use part of each. I don't think this can be avoided while using pre-determined list lengths.                                
    kipping_short = spst.beta(0.697, 3.27).rvs(size=int(num_points))
    kipping_long = spst.beta(1.12, 3.09).rvs(size=int(num_points))

    
    for i in range(num_points):
    
        per = per_list[i]
        
        if per <= 382.3:
            e = kipping_short[i]
            
        elif per >= 382.3:
            e = kipping_long[i]

        
        if e > 0.99:
            e = 0.99
        e_list[i] = e
    
    return e_list
    
    

def contour_levels(prob_array, sig_list, t_num = 1e3):
    """
    Contour drawing method taken from 
    https://stackoverflow.com/questions/37890550/python-plotting-percentile-contour-lines-of-a-probability-distribution
    This function takes a 2-D array of probabilities and returns a 1-D array 
    of the probability values p1, p2, and p3 such that a contour drawn through 
    the probabilities equal to p1 would encompass 68% of the total probability. 
    Similarly, a contour through the probabilities equaling p2 would encompass
    95% of the total probability, and 99.7% for p3. The user can specify a subset
    of [p1, p2, p3] with the sig_list argument.
    
    Arguments:
        prob_array (2D array of floats): Normalized probability array
        sig_list (list of ints): Any combination of [1,2,3] to indicate
                                 the desired probability encompassed by
                                 each contour
    
    Returns:
        t_contours (list of floats): List giving the probability value at which
                                     to draw contours. Will have the same length
                                     as sig_list
    """


    # An array of probabilites from 0 to prob_max in rate_array
    t = np.linspace(0, np.array(prob_array).max(), int(t_num))

    # (prob_array >= t[:, None, None]) is a 3D array of shape (array_num, array_num, t_num). Each (array_num, array_num) layer is a 2D array of bool values indicating which values in prob_array are greater than the value of the given t step.
    # Multiplying this 3D array of bools by prob_array replaces the bools with the array value if the bool is T and 0 if the bool is F.
    # Finally, sum along the array_num axes to get a single list of values, each with the total summed probability in its array.
    # integral is a 1D array of floats. The ith float is the sum of all probabilities in prob_array greater than the ith probability in t. Because prob_array is normalized, integral starts with 1 (all values are greater than 0, which is the first value in t) and ends with 0 (no values are greater than the greatest value).
    # Example t: t = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6] (0.6 is the maximum array value)
    # Example integral: integral = [1, 0.95, 0.8, 0.66, 0.43, 0.1, 0]

    integral = ((prob_array > t[:, None, None])*prob_array).sum(axis=(1,2))

    # Now create a function that takes integral as the x (not the y) and then returns the corresponding prob value from the t array. Interpolating between integral values allows me to choose any enclosed total prob. value (ie, integral value) and get the corresponding prob. value to use as my contour.
    # Use zero-order spline to address the interpolation issues that come with highly concentrated probability regions.
    f = sp.interpolate.interp1d(integral, t, kind='zero')

    contour_list = []
    prob_list = [0.68, 0.95, 0.997]
    
    for i in sig_list:
        contour_list.append(prob_list[i-1])

    # The plt.contourf function requires at least 2 contour values. So if we want just one contour, include another contour that encompasses almost exactly the same total probability.
    if len(sig_list) == 1:
        contour_list.append(contour_list[0]-1e-4)

    # Make sure contour_list is in descending order
    t_contours = f(np.array(sorted(contour_list, reverse=True)))
    
    # Make sure that no two probability contours have identical values. This generally only occurs for the imaging posterior in the approximate case, which is designed so that all pixels have either p=0 or p=some single value.
    for i in range(len(t_contours)):
      if i == 0:
          continue
      if t_contours[i] == t_contours[i-1]:
          t_contours[i] = t_contours[i-1]*1.001
          
    # Return t_countours, which looks like eg [0.0004, 0.0015, 0.0062]. It will be passed to matplotlib's contourf() function.
    return t_contours


def bounds_1D(prob_array, value_spaces, sig):
    """
    Given a 2D probability array, this function collapses the array along each 
    axis to find the desired confidence interval.

    Arguments:
        prob_array (2D array of floats): Square array of probabilities
        value_spaces (list of 2 2-tuples): Mass and semi-major axis limits, in the form 
                                         [(min_value_m, max_value_m), (min_value_a, max_value_a)]
        sig (int): Standard deviation limits to compute. If sig==1, compute the 68% ~ 1σ limits.

    Returns:                                       
        bounds_list (list of 2 lists): a/m values corresponding to the indices below.
        inds_sig_list (list of 2 lists): The indices on the horizontal axis where the CDF of 
                                         the collapsed prob_array reaches the upper/lower limit
                                         determined by sig.
    """
    lvls_sig_list = []
    inds_sig_list = []
    bounds_list = []

    # First compute bounds for a by collapsing along the m (aka i=0) axis
    for i in range(2):

        array_1D = prob_array.sum(axis=i)
        grid_num = len(array_1D)    

        inds_sig = CDF_indices(array_1D, [sig])[0]

        ### value_bounds is a tuple of actual values, not indices.
        ### Reverse the order of value_spaces because if we are interested in the a bounds, we collapse along m
        value_bounds = index2value(inds_sig, (0, grid_num), value_spaces[::-1][i])

        inds_sig_list.append(inds_sig)

        bounds_list.append(value_bounds)

    return bounds_list, inds_sig_list


def CDF_indices(prob_list, sig_list):
    """
    This function is the 1D analog of contour_levels(). Given a list of 
    probabilities representing a probability density function (or
    probability mass function because it is discrete), it determines 
    the interval containing a specified fraction of the total probability.
    
    Arguments:
        prob_list (array of floats): Normalized probability density function
        sig_list (list of ints): Any combination of [1,2,3] to indicate
                                 the desired probability encompassed by
                                 each set of indices. [1,2,3]-sigma
                                 correspond to [0.68, 0.95, 0.997].
    
    Returns:
        nsig_inds (list of tuples of floats): Then nth tuple in nsig_inds 
                                              gives the lower and upper
                                              indices within which is contained
                                              the probability given by the nth
                                              sigma value in sig_list.
    """
    
    # First make a list of indices in prob_list. This list is one element longer than prob_list (see below).
    ind = np.linspace(0, len(prob_list), len(prob_list)+1)

    # The input prob_list is the PDF. Use cumsum to calculate the CDF.
    CDF = np.cumsum(prob_list)
    # Insert 0 at the beginning of the cumulative sum (now the length matches ind).
    CDF = np.insert(CDF,0,0)
    # Eg: ind = [0,1,2,3,4] ; CDF = [0, 0.15, 0.4, 0.9, 1.0]
    # Matching this up with ind, we are saying that at the 0th index, we have 0 prob. At the 1st index (and after adding the 0th), we have the prob corresponding to the 1st probability sum, and so on. Depending on where the index is actually placed (I believe it's in the center of each grid block), this could incur a ~pixel-level error.
    
    # Now we have a list of indices, running from eg low mass to high mass, AND the cumulative sum at each index.
    # I want to be able to put in a cumulative probability and get out the index where the CDF attains that value.
    f = sp.interpolate.interp1d(CDF, ind)
    
    # n_sig_inds will be a list of 2-tuples. Each 2-tuple contains the indices marking the nth-sigma interval.
    # Eg, the first element might be (38.3, 65.9), which are the indices which encompass 68% of the total probability.
    nsig_inds = []
    nsig_prob_list = [0.68, 0.95, 0.997]
    
    for i in sig_list:
        prob = nsig_prob_list[i-1]
        
        # This method demands that there be equal probability excluded on both sides of the interval, which can give misleading results for multimodal or extended distributions.
        prob1 = (1-prob)/2 # Eg, (1-0.95)/2 gives the 2.5% as the first prob...
        prob2 = 1-prob1 # ... and 97.5% as the second
        
        bounds_indices = f(prob1), f(prob2)
        
        nsig_inds.append(bounds_indices)
    
    return nsig_inds
    

def index2value(index, index_space, value_space):
    """
    Converts indices to physical values.
    The axis values for a plotted array are just the array indices.
    The objective is to convert these to M and a values on a log
    scale. This function takes an index from a linear index range,
    and converts it to a parameter value in log space. 
    
    Arguments:
        index (float): Value in index space. Need not be an integer.
        index_space (tuple of floats): Index range in the form 
                                       (min_value, max_value)
        value_space (tuple of floats): Mass or semi-major axis range
                                       in the form 
                                       (min_value, max_value)
                                       
    Returns:
        value (float): Physical value corresponding to the given index
    """
    # The log base doesn't matter (it cancels) as long as it's consistent.
    index = np.array(index)

    min_index, max_index = index_space[0],  index_space[1]
    min_value, max_value = value_space[0], value_space[1]

    index_range = max_index - min_index
    log_value_range = np.log(max_value) - np.log(min_value)

    # Convert from a linear space of indices to a linear space of log(values).
    log_value = (index-min_index)*(log_value_range/index_range) + np.log(min_value)

    value = np.around(math_e**(log_value), 2) # Round to 2 decimal places

    return value


def value2index(value, index_space, value_space):
    """
    The inverse of index2value: take a value on a log scale
    and convert it to an index. 
    
    Arguments:
        value (float): Physical mass/semi-major axis value
        index_space (tuple of floats): Index range in the form 
                                       (min_value, max_value)
        value_space (tuple of floats): Mass or semi-major axis range
                                       in the form 
                                       (min_value, max_value)
    
    Returns:
        index (float): Value in index space corresponding to the
                       given physical value
    """
    min_index, max_index = index_space[0],  index_space[1]
    min_value, max_value = value_space[0], value_space[1]

    index_range = max_index - min_index
    log_value_range = np.log(max_value) - np.log(min_value)

    value_arr = np.array(value)

    index = (np.log(value_arr)-np.log(min_value))\
            *(index_range/log_value_range) + min_index

    return index
    
def min_a_and_m(trend, curv, rv_baseline, min_per, m_star):
    """
    Estimate a lower bound on the companion mass and semi-major
    axis, given the amount of RV variation observed so far.
    
    Arguments:
        trend (float, m/s/day): Linear RV trend
        curv (float, m/s/day/day): Quadratic RV curvature
        rv_baseline (float, days): Time interval over which the 
                                   RVs from which trend and curv 
                                   are calculated were taken
        min_per (float, days): The minimum period the companion
                               could have. Usually some factor larger
                               than rv_baseline, or else we would see
                               periodicity.
        m_star (float, M_Jup): Mass of host star
    
    Returns:
        min_a (float, au): Estimated minimum semi-major axis, used
                              as lower bound for sampling model sma
        min_m (float, M_Jup): Estimated minimum companion mass, used
                              as lower bound for sampling model masses
    """
    
    # Start with the minimum RV semi-amplitude K.
    # The lowest K could be is 1/2 of the current observed RV variation
    min_K = abs(0.5*(trend*rv_baseline + curv*rv_baseline**2))
    
    # Now calculate Msini with minimum period and K amplitude
    # Make sure m_star is in solar masses
    # e=0 for simplicity, though mass could be lower if e were very high
    # High-e companion near min_per is an edge case that can be handled separately
    min_m = rv.utils.Msini(min_K, min_per, m_star*M_jup/M_sun, 
                           0, Msini_units='jupiter')
                           
    min_m = max(min_m, 0.1) # min_mass cannot be 0 or logarithmic spacing breaks
    
    min_a = rv.utils.semi_major_axis(min_per, ((m_star + min_m)*(M_jup/M_sun)))
    
    return min_a, min_m
    
    