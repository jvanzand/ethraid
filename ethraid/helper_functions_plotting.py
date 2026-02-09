import os
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.patches as ptch

from ethraid.compiled import helper_functions_general as hlp

two_pi = 6.283185307179586
G = 2.824760877012879e-07 # c.G.cgs.value*(1/c.au.cgs.value)**3 * (c.M_jup.cgs.value) * (24*3600)**2

def scatter_companion(scatter_plot, grid_num_ext, a_lim_plot, m_lim_plot):
    """
    Simple function to calculate where to plot an expected/known companion.
    
    Arguments:
        scatter_plot (tuple of floats): (semi-major axis, mass) values at
                                        which to plot a companion. Sma is
                                        in AU and mass is in M_Jup.
        grid_num_ext (int): Dimension of square plotting array, including
                           extension for ruled-out regions
        a_lim_plot (tuple of floats, au): Semi-major axis limits to plot 
                                          over, in the form (a_min, a_max)
        m_lim_plot (tuple of floats, M_jup): Mass limits as (m_min, m_max)
    
    Returns:
        sep_ind, mp_ind (tuple of floats): Indices corresponding to given
                                           separation and mass values.
    """
    
    sma, mass = scatter_plot
    
    sep_ind = hlp.value2index(sma, (0, grid_num_ext-1), a_lim_plot)
    mp_ind  = hlp.value2index(mass, (0, grid_num_ext-1), m_lim_plot)
    

    scatter_a, scatter_m = scatter_plot
    c1 = scatter_a<a_lim_plot[0]
    c2 = scatter_a>a_lim_plot[1]
    c3 = scatter_m<m_lim_plot[0]
    c4 = scatter_m>m_lim_plot[1]

    if c1 or c2 or c3 or c4:
        print("helper_functions_plotting.py: WARNING -\n"\
              "        Attempting to plot companion coordinates outside plotting bounds.\n"\
              "        Increase plotting bounds to encompass companion parameters.")
    
    return sep_ind, mp_ind
        

def period_lines(m_star, a_lim, m_lim, a_lim_plot, m_lim_plot, grid_num_ext, n, how='tot'):
    """
    Calculates indices of lines of constant period to plot over 2-D posterior.
    
    Arguments:
        m_star (float, M_jup): Stellar mass
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max)
        m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
        a_lim_plot (tuple of floats, au): Semi-major axis limits including
                                          buffer for ruled-out regions. 
                                          The interval a_lim is contained 
                                          in a_lim_plot.
        m_lim_plot (tuple of floats, M_jup): Companion mass limits including
                                             buffer for ruled-out regions. 
                                             The interval m_lim is contained 
                                             in m_lim_plot.
        grid_num_2D (int): Shape of square probability array that is used
                           for plotting
        n (int): How many lines of constant period to plot
        how (str): Either 'tot' or 'Gaia'. Determines which harmonics to plot
    
    Returns:
        trimmed_a_inds_list (list of floats): List of indices for plotting a values
        trimmed_m_inds_list (list of floats): List of indices for plotting m values
        fmt (str): Format of plotted line
    """
    
    a_min, a_max = a_lim
    m_min, m_max = m_lim
    
    ######## Adding lines of constant period ##########
    hip_times  = [Time(1989.85, format='decimalyear').jd, Time(1993.21, format='decimalyear').jd]       
    #https://www.cosmos.esa.int/web/hipparcos/catalogue-summary

    gaia_times = [Time('2014-07-25', format='isot').jd, Time('2017-05-28', format='isot').jd] 
    #https://www.cosmos.esa.int/web/gaia/earlydr3

    # Time between the midpoints of the two missions
    baseline_days = ((gaia_times[1] + gaia_times[0])/2 - (hip_times[1] + hip_times[0])/2)
    gaia_baseline_days = gaia_times[1] - gaia_times[0]

    # Masses on the vertical axis. Simply log-spaced.
    const_per_m_list = np.logspace(np.log10(m_min), np.log10(m_lim[1]), 50)
    const_per_m_inds = hlp.value2index(const_per_m_list, (0, grid_num_ext-1), m_lim_plot)
    
    # Space the lines at harmonics of the mission baseline. Harmonics of 1) the Gaia mission baseline only and 2) the baseline between Hipparcos and Gaia, are relevant.
    if how == 'tot':
        baseline = baseline_days
        fmt = '--k'
    elif how == 'gaia':
        baseline = gaia_baseline_days
        fmt = '--r'
    else:
        raise Exception('helper_functions_plotting.py: Argument "how" must be either "tot" or "gaia"')

    # Lines of constant period for p = baseline_days/n
    #
    trimmed_a_inds_list = []
    trimmed_m_inds_list = []
    for i in range(n):
        const_per_a_list = constant_per_a(const_per_m_list, baseline/(i+1), m_star)
        const_per_a_inds = hlp.value2index(const_per_a_list, (0, grid_num_ext-1), a_lim_plot)
        
        values_in_bounds = np.where((a_lim[0] < const_per_a_list)&(const_per_a_list < a_lim[1]))
        
        
        trimmed_a_inds_list.append(const_per_a_inds[values_in_bounds])
        trimmed_m_inds_list.append(const_per_m_inds[values_in_bounds])

    
    return trimmed_a_inds_list, trimmed_m_inds_list, fmt

def constant_per_a(m, per, m_star):
    """
    Function to help draw lines of constant period on 
    the final plot. Rearranges Kepler's 3rd law to find 
    how semi-major axis a varies with period, companion 
    mass, and stellar mass.

    Intended usage: Calculate an array of a values for a fixed per
                    and m_star and an array of companion masses.
            
    Arguments:
            m (list of floats, M_J): companion masses
            per (float, days): Fixed orbital period
            m_star (float, M_J): stellar mass

    Returns:
            a (list of floats, au): Semi-major axis values (au)
                                    which produce orbits with
                                    period = per, given the
                                    masses m.
    """

    a = ((per/two_pi)**2*G*(m+m_star))**0.3333333333333


    return a
        


def sma2angsep(sma, d_star):
    """
    Converts separation in AU into
    separation in arcsec
    
    Arguments:
        sma (float, AU): Projected physical separation of
                         a companion from its host
        d_star (float, AU): Distance of host star from Earth
    
    Returns:
        ang_sep (float, arcsec): Projected angular host-
                                 companion separation
    """
    pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
    d_star_pc = np.array(d_star)/np.array(pc_in_au)
    ang_sep = sma/d_star_pc
    
    return ang_sep

def angsep2sma(ang_sep, d_star):
    """
    Inverse of sma2angsep.
    Converts separation in AU into
    separation in arcsec
    
    Arguments:
        ang_sep (float, arcsec): Projected angular host-
                                 companion separation
        d_star (float, AU): Distance of host star from Earth
    
    Returns:
        sma (float, AU): Projected physical separation of
                         a companion from its host
    """
    pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
    d_star_pc = np.array(d_star)/np.array(pc_in_au)
    sma = ang_sep*d_star_pc
    
    return sma   
        
        
        