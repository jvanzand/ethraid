import sys
import numpy as np
import scipy as sp
import scipy.stats as spst

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

def make_arrays(double m_star, tuple a_lim, tuple m_lim, int grid_num, int num_points, str e_dist):
    """
    Create the parameter arrays which will be used for the RV and astrometry posteriors.
    
    Arguments:
        m_star (float, M_sun): Mass of host star
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max)
        m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
        grid_num (int): Dimensions of square (a,m) grid
        num_points (int): Number of random orbital models to simulate
        e_dist (str): Which eccentricity distribution to draw from.
                        See ecc_dist() function below.
    
    Returns:
        a_list, m_list, per_list, e_list, 
        i_list, om_list, M_anom_0_list (1D numpy arrays, len = num_points):
                                        Lists of randomly sampled semi-major axis, mass,
                                        eccentricity, inclination, argument of
                                        periastron, and initial mean anomaly. Do not
                                        sample directly in period.
        a_inds, m_inds (1D numpy arrays of ints, len = num_points): Grid position where each 
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
                                    om_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    a_bins = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    m_bins = np.ndarray(shape=(num_points,), dtype=np.float64)

    cdef long [:] a_inds, m_inds

    tp = 0
    a_min = a_lim[0]
    a_max = a_lim[1]
    m_min = m_lim[0]
    m_max = m_lim[1]

    # These are the "full" semi-major axes of the orbit, ie the sma of the ellipse traced by the 1-body solution to the 2-body problem. a = a_planet+a_star
    a_list = spst.loguniform.rvs(a_min, a_max, size=num_points)
    m_list = spst.loguniform.rvs(m_min, m_max, size=num_points)
    
    a_prior = spst.loguniform.pdf(a_list, a_min, a_max) # Find the PDF value of each value in a_list
    m_prior = spst.loguniform.pdf(m_list, m_min, m_max) # Find the PDF value of each value in m_list

    # Match up a_list and m_list and get the period for each pair (in days).
    # Calculate this now to avoid having to do it twice for RVs and astrometry.
    per_list = P_list(a_list, m_list, m_star)
    
    # Eccentricities drawn from specified distribution.
    e_list, e_prior = ecc_dist(m_list, per_list, num_points, dist=e_dist)

    cosi_list = np.random.uniform(0, 1, num_points)
    i_list = np.arccos(cosi_list)
    i_prior = np.sin(i_list) # PDF of inclination is prop to sin(i)

    # Mean anomaly, uniformly distributed. This represents M at the beginning of the Hipparcos epoch for BOTH RVs and astrometry. Use this to solve for True anomaly.
    M_anom_0_list = np.random.uniform(0, two_pi, num_points)

    # Arguments of peri of the companion, uniformly distributed
    om_list = np.random.uniform(0, two_pi, num_points)

    # Breaking up the (a, M) parameter space into grid_num x grid_num.
    # Note: grid_num+1 bin dividers means grid_num+2 bins (including below/above min/max), but the first/last bins are never used because we sample between the min and max values. So if grid_num=100, there are bins 0 up to 101. But only bins 1 to 100 have values in them.
    a_bins = np.logspace(np.log10(a_min), np.log10(a_max), grid_num+1)
    m_bins = np.logspace(np.log10(m_min), np.log10(m_max), grid_num+1)

    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)
    
    tot_prior = a_prior*m_prior*e_prior*i_prior # Multiply priors to obtain total prior (note: not log)

    return a_list, m_list, per_list, e_list, i_list,\
           om_list, M_anom_0_list, a_inds, m_inds, tot_prior

def tot_list(double [:] rv_post_list, double [:] astro_post_list, 
             double [:] imag_post_list, int num_points):
    """
    Start with 3 1D lists and add them element-wise.
    
    Arguments:
        rv_post_list, astro_post_list, 
        imag_post_list (1D numpy arrays, len=num_points):
            Lists of log-likelihoods of the rv/astrometry/imaging data
            conditioned on a given model. Then ith log-likelihood in each
            list corresponds to the same model.
        num_points (int): Length of above lists. Also equal to the
                          total number of models run.
    
    Returns:
        tot_list (1D numpy array, len=num_points):
            The log-likelihoods of all data sets, conditioned on a given
            model. The ith element is the sum of the ith elements of the 
            three above lists.
    """
    cdef int i
    cdef double prob
    cdef double [:] tot_list = np.zeros(shape=(num_points), dtype=np.float64)

    for i in range(num_points):

       log_prob = rv_post_list[i]+astro_post_list[i]+imag_post_list[i]

       tot_list[i] += log_prob
    
    return tot_list



def post_tot_approx_imag(double [:] tot_list, double [:,:] post_imag,
                         long [:] a_inds, long [:] m_inds, int grid_num):
            
    """
    Start with 2 1D lists and multiply them element-wise, THEN form 
    the result into a 2D array. This function is for the total posterior;
    the individual RV, astrometry, and imaging posteriors are handled by the 
    post_single() function below.
    
    Arguments:
        tot_list (np array of floats, len=num_points): List of RV/astro data log-likelihoods
        post_imag (2D np array of floats, dim=num_points x num_points): Array of imaging data log-likelihoods
        a_inds, m_inds (numpy arrays of ints, len = num_points): Grid position where each 
                                        (a, m, per, e, i, om, M_anom_0) model will 
                                        be placed, based on the model's 
                                        a and m values
        grid_num (int): Dimension of square (a,m) grid
                                        
    Returns:
        tot_prob_array (numpy array, dim = grid_num x grid_dum): 2-D array of binned
                       (aka marginalized) posterior probabilities. Note, the binning
                       process itself is what applies the priors, converting the
                       individual likelihoods into posterior probabilities.
    """

    cdef int i, a_i, m_i
    cdef double [:,:] rv_ast_array = np.zeros(shape=(grid_num,grid_num), dtype=np.float64)
    
    rv_ast_array = post_single(tot_list, a_inds, m_inds, grid_num)
    
    
    # Not cdef because then I can't change it from memview to np array
    tot_prob_array = np.zeros((grid_num, grid_num))
    
    # post_imag is not a list like RVs and astrometry above. It is input to this function as a 2D array. This is because it's a lot easier to calculate. It would take way longer if we calculated a length-1e6 (or 1e8) list instead of a 100x100 array.
    for i in range(grid_num):
        for j in range(grid_num):
            tot_prob_array[i,j] += rv_ast_array[i,j]*post_imag[i,j]
    

    return np.array(tot_prob_array)/np.array(tot_prob_array).sum()


def post_single(double [:] log_lik_list, long [:] a_inds, long [:] m_inds, int grid_num):
    """
    Form a list of log-likelihoods into a 2D array of (non-log) posterior probabilities.
    
    Arguments:
        log_lik_list (list/array): List of log-likelihoods to be be reshaped
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

    size = log_lik_list.shape[0]
    for i in range(size):
        # The lowest index in a_inds or m_inds is 1 because the lower limit of the sampling range (a_min or m_min) is also the smallest bin divider. So no sampled values are small enough to get an index of 0. Still, we want the bin between a_min and the next divider to be the 0th bin, so subtract 1 off of each bin value.
        a_i = a_inds[i]-1
        m_i = m_inds[i]-1

        prob_array[m_i, a_i] += np.exp(log_lik_list[i])
        
    
    return np.array(prob_array)/np.array(prob_array).sum()


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
    semimajor axis, planet mass, and the host star mass.
    
    Arguments:
        a (float, au): Semi-major axis
        m_planet (float, M_jup): Companion mass
        m_star (float, M_jup): Stellar mass
    
    Returns:
        per (float, days): Orbital period
    """
    
    cdef double per
    
    per = two_pi * a**(1.5) * (G*(m_planet+m_star))**(-0.5)
    
    return per


def ecc_dist(double [:] m_list, double [:] per_list, int num_points, dist='piecewise'):
    """
    Sample a random eccentricity whose distribution is based on a and m.
    Kipping(2013) advocates two distributions for P below and above 382.3 days.
    
    Arguments:
        m (list of floats, M_jup): List of companion masses to choose
                                   ecc sub-distribution
        per_list (list of floats, days): List of companion periods
        num_points (int): Number of random eccentricities to generate
        dist (str): Desired eccentricity distribution
                    Options: ['zero', 'uniform', 'kipping', 'piecewise']
    
    Returns:
        e_list (list of floats): List of randomly sampled eccentricities
    """
    cdef int i
    cdef double m, per, e
    
    cdef np.ndarray[double, ndim=1] e_list = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_prior = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_sample1 = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_sample2 = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_sample3 = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_sample4 = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_prior1 = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_prior2 = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_prior3 = np.ndarray(shape=(int(num_points)), dtype=np.float64),\
                                    e_prior4 = np.ndarray(shape=(int(num_points)), dtype=np.float64)
    
    # Fix e=0
    if dist=='zero':
        e_list = np.zeros(num_points)
        e_prior = np.ones(int(num_points))
        
    # Draw e uniformly between 0 and 0.99
    elif dist=='uniform':
        e_list = spst.uniform(0, 0.99).rvs(size=int(num_points))
        e_prior = np.ones(int(num_points))
        
    # Use Kipping (2013) distribution for all companions
    elif dist=='kipping':
        # Note that I make lists that are each num_points long, so I only use part of each. I don't think this can be avoided while using pre-determined list lengths.                                
        e_sample1 = spst.beta(0.697, 3.27).rvs(size=int(num_points))
        e_sample2 = spst.beta(1.12, 3.09).rvs(size=int(num_points))
        
        e_prior1 = spst.beta(0.697, 3.27).pdf(e_sample1)
        e_prior2 = spst.beta(1.12, 3.09).pdf(e_sample2)

    
        for i in range(num_points):
            per = per_list[i]
        
            if per <= 382.3:
                e = e_sample1[i]
                e_pri = e_prior1[i]
            
            elif per > 382.3:
                e = e_sample2[i]
                e_pri = e_prior2[i]

            if e > 0.99:
                e = 0.99
            e_list[i] = e
            e_prior[i] = e_pri
    
    # Use Kipping for planetary masses, Bowler (2020) for BDs, and Raghavan (2010) for stars
    # Why not determine based on *separation* instead of mass? Preliminary CLS analysis suggests e depends more on mass than on a.
    elif dist=='piecewise':
        e_sample1 = spst.beta(0.697, 3.27).rvs(size=int(num_points))
        e_sample2 = spst.beta(1.12, 3.09).rvs(size=int(num_points))
        e_sample3 = spst.beta(2.30, 1.65).rvs(size=int(num_points))
        e_sample4 = spst.uniform(0.1, 0.7).rvs(size=int(num_points)) # (0.1, 0.7) --> sample from [0.1,0.8]
        
        e_prior1 = spst.beta(0.697, 3.27).pdf(e_sample1)
        e_prior2 = spst.beta(1.12, 3.09).pdf(e_sample2)
        e_prior3 = spst.beta(2.30, 1.65).pdf(e_sample3)
        e_prior4 = spst.uniform(0.1, 0.7).pdf(e_sample4) # (0.1, 0.7) --> sample from [0.1,0.8]
        
        for i in range(num_points):
            m = m_list[i]
            
            if m <= 13:
                per = per_list[i]
                
                if per <= 382.3:
                    e = e_sample1[i]
                    e_pri = e_prior1[i]
            
                elif per > 382.3:
                    e = e_sample2[i]
                    e_pri = e_prior2[i]
                    
            elif m > 13 and m <= 80:
                e = e_sample3[i]
                e_pri = e_prior3[i]
            
            elif m > 80:
                e = e_sample4[i]
                e_pri = e_prior4[i]
            
            if e > 0.99:
                e = 0.99
            e_list[i] = e
            e_prior[i] = e_pri
            
    else:
        raise Exception('Error: e_dist must be one of the options supported by the ecc_dist function.')
    
    
    return e_list, e_prior
    
    

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
    Given a normalized 2D probability array, this function collapses the array along each 
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
    This function is the 1D analog of contour_levels(). Given a normalized
    list of probabilities representing a probability density function (or
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
    # Eg: PDF = [0.15, 0.25, 0.5, 0.1] ; ind = [0,1,2,3,4] ; CDF = [0, 0.15, 0.4, 0.9, 1.0]
    # Matching this up with ind, we are saying that at the 0th index, we have 0 prob. At the 1st index (and after adding the 0th), we have the prob corresponding to the 1st probability sum, and so on.
    
    # Now we have a list of indices, running from eg low mass to high mass, AND the cumulative sum at each index.
    # Goal: to be able to put in a cumulative probability and get out the index where the CDF attains that value.
    f = sp.interpolate.interp1d(CDF, ind)
    
    # n_sig_inds will be a list of 2-tuples. Each 2-tuple contains the indices marking the nth-sigma interval.
    # Eg, the first element might be (38.3, 65.9), which would be the indices which encompass 68% of the total probability.
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
    
    try:
        import radvel as rv
    except ModuleNotFoundError as err:
        raise Exception('helper_functions_general.min_a_and_m: Radvel 1.3.8 must be installed to use this function.')
    
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
    
    