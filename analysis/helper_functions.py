import numpy as np
import scipy as sp
import astropy.constants as c
import astropy.units as u

import radvel as rv
from c_kepler import _kepler as ck

# import matplotlib.pyplot as plt
        
        
def P(a, Mtotal):
    """
    Uses Kepler's third law to find the period of a planet (in days) given its 
    semimajor axis and the total mass of the system.
    
    a (au): semi-major axis
    Mtotal (Msun): Mass of star + mass of object
    """
    
    Mtotal_grams = Mtotal*c.M_sun.cgs.value
    a_cm = a*c.au.cgs.value
    
    P_sec = np.sqrt((2*np.pi)**2*a_cm**3/(c.G.cgs*Mtotal_grams).value)
    
    P_days = P_sec/(24*3600)
    
    return P_days

def rot_matrix(i, om, Om):
    """
    This is P3*P2*P1 from Murray & Dermott. It is not given explicitly in the text. They multiply it immediately by r*[cos(f), sin(f), 0]
    because this gives the projection of position onto the sky. However, we also need the projection of velocity, so we need the matrix
    pre-multiplication by the position vector.
    """
    row_1 = [np.cos(Om)*np.cos(om) - np.sin(Om)*np.cos(i)*np.sin(om), 
            -np.sin(om)*np.cos(Om) - np.sin(Om)*np.cos(i)*np.cos(om), 
            np.sin(Om)*np.sin(i)]
            
    row_2 = [np.sin(Om)*np.cos(om) + np.cos(Om)*np.cos(i)*np.sin(om), 
            -np.sin(om)*np.sin(Om) + np.cos(Om)*np.cos(i)*np.cos(om), 
            -np.cos(Om)*np.sin(i)]
            
    row_3 = [np.sin(i)*np.sin(om), np.sin(i)*np.cos(om), np.cos(i)]
    
    rot_matrix = np.array([row_1, row_2, row_3])
    
    return rot_matrix
    

def r(nu, a, e):
    """
    
    Equation of an ellipse (Murray & Dermott equation 2.20).
    Arguments:
    
        nu (radians): True anomaly
        a (distance): Semi-major axis of ellipse. Choice of a determines what output r represents. For example, if a is the semi-major
                                       axis of one planet's orbit, then r represents that planet's distance from barycenter as a function of nu. On the
                                       other hand, if a is the SA of the test mass μ's orbit, then r is r1+r2 as a function of nu, where r1 (r2) is the 
                                       distance of m1 (m2) from the system barycenter in the 2-body (m1 & m2) frame.
        e (unitless): Eccentricity
    
    returns:
        r (same as a): Distance of particle from barycenter along its orbit
    """
    num = a*(1-e**2)
    denom = 1 + e*np.cos(nu)
    
    return num/denom
    
def v_vec(a, per, e, nu):
    """
    Uses Murray & Dermott equation 2.36. r_dot is not what we want because it doesn't capture the velocity perpendicular to the radial vector.
    Instead, v is the total velocity of the object. M&D doesn't actually give v vector explicitly, but I believe it's v_vec = [x_dot, y_dot, 0]. 
    
    All of the input arrays must have compatible shapes.
    """
    n = 2*np.pi/per
    
    x_dot = -n*a / np.sqrt(1-e**2) * np.sin(nu)
    y_dot = +n*a / np.sqrt(1-e**2) * (e + np.cos(nu))
    
    # To get the proper shape vector at the end, we need our 0 element to be an array with matching shape
    zero_shape = np.shape(y_dot)
    
    v_vec = np.array([x_dot, y_dot, np.zeros(zero_shape)])
    
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
    
def deriv(y_array, x_array):
    """
    Computes an array of first derivative values. y_array and x_array must be the same length.
    The derivative at each point is estimated by finding the slope between the point after and 
    the point before. This means that the array of derivatives will have 2 fewer elements than
    y_array and x_array.
    """
    
    y_array = np.array(y_array)
    x_array = np.array(x_array)
    
    deriv_array = (y_array[2:] - y_array[0:-2])/(x_array[2:] - x_array[0:-2])
    
    return deriv_array 

def gamma(a, Mp, per, e, i, om, M_anom):
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
    
    Mp = Mp * c.M_jup.cgs.value
    G = c.G.cgs.value
    a = a*c.au.cgs.value
    
    # E = kepler(M_anom, np.array([e for i in range(len(M_anom))]))
    # M_anom = M_anom[:, None] * np.linspace(0, 100, int(1e2))
    # e = e[:, None]
    
    
    # E = kepler(M_anom, e)
    print('using c-solver')
    E = ck.kepler_array(M_anom, e)
    
    # E = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(nu/2))
    
    nu = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))
    
    # Differentiate Kepler's equation in time to get E_dot
    # Note that E_dot has units of (1/per), where [per] is days. Therefore [gamma_ddot] = m/s/d^2 
    E_dot = (2*np.pi/per)/(1-e*np.cos(E))
    nu_dot = (1+np.tan(nu/2)**2)**-1 * np.sqrt((1+e)/(1-e)) * np.cos(E/2)**-2 * E_dot
    
    # Convert prefac units from cm/s^2 to m/s/day
    # Negative just depends on choice of reference direction. I  am being consistent with radvel rv_drive function.
    prefac = -(Mp*G*np.sin(i))/(a**2*(1-e)) * (1/100) * (24*3600)
    

    
    gd_t1 = (1+np.cos(nu))/(1+np.cos(E))
    gd_t2 = np.sin(nu+om)/(1-e*np.cos(E))

    
    gamma_dot = prefac*gd_t1*gd_t2
    
    gd_t1_dot = ((1+np.cos(nu))*np.sin(E) * E_dot - (1+np.cos(E))*np.sin(nu)*nu_dot) / (1+np.cos(E))**2
    gd_t2_dot = ((1-e*np.cos(E))*np.cos(nu+om) * nu_dot - np.sin(nu+om)*e*np.sin(E)*E_dot) / (1-e*np.cos(E))**2
    
    # gdd_t1 = gd_t2 * ((1+np.cos(nu))*np.sin(E) * E_dot - (1+np.cos(E))*np.sin(nu)*nu_dot) / (1+np.cos(E))**2
    # gdd_t2 = gd_t1  * ((1-e*np.cos(E))*np.cos(nu+om) * nu_dot - np.sin(nu+om)*e*np.sin(E)*E_dot) / (1-e*np.cos(E))**2
    
    gdd_t1 = gd_t2 * gd_t1_dot
    gdd_t2 = gd_t1 * gd_t2_dot
    
    gamma_ddot = prefac*(gdd_t1+gdd_t2)

    return gamma_dot, gamma_ddot

def gamma_T(a, Mp, per, e, i, om, nu):
    """
    ## This function does NOT solve Kepler's equation! It requires nu to be provided directly.
    Function to analytically calculate the first and second derivatives of the RV curve at a given point in the orbit.
    Mp is expected in Jupiter masses.
    a is expected in au
    per is expected in days
    
    Returns:
    gamma_dot (m/s/d)
    gamma_ddot (m/s/d^2)
    """
    
    Mp = Mp * c.M_jup.cgs.value
    G = c.G.cgs.value
    a = a*c.au.cgs.value
    
    E = 2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(nu/2))
    
    E_dot = (2*np.pi/per)/(1-e*np.cos(E))
    nu_dot = (1+np.tan(nu/2)**2)**-1 * np.sqrt((1+e)/(1-e)) * np.cos(E/2)**-2 * E_dot
    
    # Convert prefac units from cm/s^2 to m/s/day
    # Negative just depends on choice of reference direction. I  am being consistent with radvel rv_drive function.
    prefac = -(Mp*G*np.sin(i))/(a**2*(1-e)) * (1/100) * (24*3600)
    

    
    gd_t1 = (1+np.cos(nu))/(1+np.cos(E))
    gd_t2 = np.sin(nu+om)/(1-e*np.cos(E))

    
    gamma_dot = prefac*gd_t1*gd_t2
    
    gd_t1_dot = ((1+np.cos(nu))*np.sin(E) * E_dot - (1+np.cos(E))*np.sin(nu)*nu_dot) / (1+np.cos(E))**2
    gd_t2_dot = ((1-e*np.cos(E))*np.cos(nu+om) * nu_dot - np.sin(nu+om)*e*np.sin(E)*E_dot) / (1-e*np.cos(E))**2
    
    # gdd_t1 = gd_t2 * ((1+np.cos(nu))*np.sin(E) * E_dot - (1+np.cos(E))*np.sin(nu)*nu_dot) / (1+np.cos(E))**2
    # gdd_t2 = gd_t1  * ((1-e*np.cos(E))*np.cos(nu+om) * nu_dot - np.sin(nu+om)*e*np.sin(E)*E_dot) / (1-e*np.cos(E))**2
    
    gdd_t1 = gd_t2 * gd_t1_dot
    gdd_t2 = gd_t1 * gd_t2_dot
    
    gamma_ddot = prefac*(gdd_t1+gdd_t2)
    
    return gamma_dot, gamma_ddot
    
def kepler(Marr, eccarr):
    """
    Solve Kepler's Equation using a modified Newton-Raphson method
    Args:
        Marr (array): input Mean anomaly
        eccarr (array): eccentricity
    Returns:
        array: eccentric anomaly
    
    From Radvel package. Adapted here to accept non-equal-length Marr and eccarr
    """
    # if not isinstance(Marr, np.ndarray):
    #     Marr = np.array([Marr])
    #     eccarr = np.array([eccarr for i in range(len(Marr))])
    #
    #
    # if isinstance(eccarr, float):
    #     eccarr =np.array([eccarr for i in range(len(Marr))])

        
    
    conv = 1.0e-12  # convergence criterion
    k = 0.85


    Earr = Marr + np.sign(np.sin(Marr)) * k * eccarr  # first guess at E

    # fiarr should go to zero when converges
    fiarr = ( Earr - eccarr * np.sin(Earr) - Marr)
    
    convd = np.where(np.abs(fiarr) > conv)[0]  # which indices have not converged
    nd = len(convd)  # number of unconverged elements
    count = 0

    while nd > 0:  # while unconverged elements exist
        count += 1
        
        M = Marr[convd]  # just the unconverged elements ...
        ecc = eccarr[convd]
        E = Earr[convd]
        
        fi = fiarr[convd]  # fi = E - e*np.sin(E)-M    ; should go to 0
        fip = 1 - ecc * np.cos(E)  # d/dE(fi) ;i.e.,  fi^(prime)
        fipp = ecc * np.sin(E)  # d/dE(d/dE(fi)) ;i.e.,  fi^(\prime\prime)
        fippp = 1 - fip  # d/dE(d/dE(d/dE(fi))) ;i.e.,  fi^(\prime\prime\prime)
        
        ##
        # fipppp = -fipp
        # fi5 = -fippp
        ##

        # first, second, and third order corrections to E
        d1 = -fi / fip
        d2 = -fi / (fip + d1 * fipp / 2.0)
        d3 = -fi / (fip + d2 * fipp / 2.0 + d2 * d2 * fippp / 6.0)
        ##
        # d4 = -fi / (fip + d3 * fipp / 2.0 + d3 * d3 * fippp / 6.0 + d3*d3*d3*fipppp / 24.0)
        # d5 = -fi / (fip + d4 * fipp / 2.0 + d4 * d4 * fippp / 6.0 + d4*d4*d4*fipppp / 24.0 + d4*d4*d4*d4*fi5 / 120.0)
        ##
        E = E + d3
        Earr[convd] = E
        fiarr = ( Earr - eccarr * np.sin( Earr ) - Marr) # how well did we do?
        convd = np.abs(fiarr) > conv  # test for convergence
        nd = np.sum(convd is True)

    if Earr.size > 1:
        return Earr
    else:
        return Earr[0]
        
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
    t = np.linspace(0, prob_array.max(), int(t_num))

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
        
        # This gives only the 2-sigma levels, so that we get the 2-sigma limits at the end
        sig2 = contour_levels_1D(array_1D, [2])[0]
        
        # Interpolate between the points to get a finer spacing of points. This allows for more precise parameter estimation.
        func = sp.interpolate.interp1d(range(grid_num), array_1D)
        
        # Array over the same interval, but spaced (probably) more finely
        fine_array = np.linspace(0, grid_num-1, int(interp_num))
        
        # This is analogous to the original array_1D, but finer
        interp_vals = func(fine_array)
        
        
        #### For debugging purposes
        # plt.plot(fine_array, interp_vals)
        # plt.show()
        # print(i, np.where(abs(interp_vals - sig2) < 1e-2*sig2)[0])
        ####
        
        # This is a shaky step. I'm just looking for places where the function value is really close to the probability corresponding to 2-sigma. But from what I can tell, this will fall apart for multimodal distributions, and maybe in other cases too. I use the 'take' method to pick out the first and last indices.

        inds_2sig = np.where(abs(interp_vals - sig2) < 1e-2*sig2)[0].take((0,-1))
        
        # value_bounds is a tuple of actual values, not indices
        value_bounds = index2value(inds_2sig, (0, interp_num-1), value_spaces[::-1][i])
        
        bounds_list.append(value_bounds)
        
    return bounds_list


def value2index(value, index_space, value_space):
    """
    The inverse of index2value: take a value on a log scale and convert it to an index. index_space and value_space are expected
    as tuples of the form (min_value, max_value).
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
    The axis values for a plotted array are just the array indices. I want to convert these to Msini and a values, and on a log
    scale. This function takes a single index from a linear index range, and converts it to a parameter value in log space.
    index_space and value_space are expected as tuples of the form (min_value, max_value). index is in the range of index_space.
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


def multiple_period_plotter(a_list, m_list, t_len = 200):
    
    fig, axs = plt.subplots(2, len(a_list), squeeze=False)
    fig.tight_layout(h_pad = 4, w_pad = 4)
    
    y_titles = [r'$\dot\gamma$ $(m/s^2)$', r'$\ddot\gamma$ $(m/s^3)$']
    fontsize = 14
    
    for i in range(len(a_list)):
        
        a = a_list[i]
        m = m_list[i]
        
        per = P(a, (m_star+m*(c.M_jup/c.M_sun).value))
        k = rv.utils.semi_amplitude(m, per, (m_star+m*(c.M_jup/c.M_sun).value), e, Msini_units='jupiter')

        orbel = [per, tp, e, om, k]


        t_list = np.linspace(0, per, t_len) # List of times spanning one full orbit. Dates are arbitrary so I start at jd=0 for ease

        rv_list = rv.kepler.rv_drive(t_list, orbel)

        gamma_dot_array = deriv(rv_list, t_list)[1:-1]
        gamma_dotdot_array = deriv(deriv(rv_list, t_list), t_list[1:-1])
        
        axs[0, i].plot(t_list[2:-2], gamma_dot_array)
        axs[0, i].plot(t_list, [gamma_dot + gamma_dot_err for i in t_list], c='r')
        axs[0, i].plot(t_list, [gamma_dot - gamma_dot_err for i in t_list], c='r')
        
        axs[1, i].plot(t_list[2:-2], gamma_dotdot_array)
        axs[1, i].plot(t_list, [gamma_dotdot + gamma_dotdot_err for i in t_list], c='r')
        axs[1, i].plot(t_list, [gamma_dotdot - gamma_dotdot_err for i in t_list], c='r')
        
        axs[0, i].set_title(r'a = {} au, Msin(i) = {} $M_{}$, $\dot\gamma$'.format(a, m, '{Jup}'), fontsize=fontsize)
        axs[0, i].set_xlabel('Time (days)', fontsize=fontsize)
        axs[0, i].set_ylabel(y_titles[0], fontsize=fontsize)
        axs[0, i].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        
        
        axs[1, i].set_title(r'a = {} au, Msin(i) = {} $M_{}$, $\ddot\gamma$'.format(a, m, '{Jup}'), fontsize=fontsize)
        axs[1, i].set_xlabel('Time (days)', fontsize=fontsize)
        axs[1, i].set_ylabel(y_titles[1], fontsize=fontsize)
        axs[1, i].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
    plt.show()
    
    
    
    
    
    
    
    
    
    
    
    