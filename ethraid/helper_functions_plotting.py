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
    # a_min, m_min = a_lim[0], m_lim[0]
    # a_max, m_max = a_lim[1], m_lim[1]
    
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
        
def marginalized_1d(star_name, post_tot, twosig_inds, a_lim, m_lim, 
                    tick_labels_a, tick_labels_m, outdir=''):
    """
    Plots and saves 2 marginalized posterior cumulative distribution function (CDF).
    The first is marginalized over mass, so it gives the semi-major axis CDF. The
    second is marginalized over semi-major axis and gives the mass CDF.
                    
    Arguments:
        star_name (str): Star name to label saved figure
        post_tot (array of floats): Posterior probability array
        twosig_inds (list of 2 lists): Each set of indices in twosig_inds
                                       encompasses 95% of the sma or mass
                                       posterior.
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max)
        m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
        tick_labels_a (list of floats, AU): Sma values to use as axis 
                                            labels
        tick_labels_m (list of floats, M_jup): Mass values to use as 
                                               axis labels
        outdir (str): Path to save plot
    
    Returns:
        None (but plots and saves 1D posteriors)
        
    """
    title_size = 30
    label_size = 25
    tick_num = 6
    tick_size = 25
    
    fig, ax = plt.subplots(1,2, figsize=(12,8))
    sma_1d = post_tot.sum(axis=0)
    mass_1d = post_tot.sum(axis=1)
    
    grid_num = np.shape(post_tot)[0]
    tick_positions_a1D = hlp.value2index(tick_labels_a, (0, grid_num-1), a_lim)
    tick_positions_m1D = hlp.value2index(tick_labels_m, (0, grid_num-1), m_lim)

    ax[0].plot(range(grid_num+1), np.insert(np.cumsum(sma_1d), 0, 0))
    plt.sca(ax[0])
    plt.xticks(tick_positions_a1D, tick_labels_a, size=tick_size)
    plt.yticks(size=tick_size)
    plt.title('Semi-major axis CDF', size=title_size)
    plt.xlabel('Companion semi-major axis (AU)', size = label_size)

    ax[0].hlines(0, 0, grid_num-1, colors='k', linestyles='solid')
    ax[0].vlines(twosig_inds[0][0], 0, 1, colors='r', linestyles='dashed')
    ax[0].vlines(twosig_inds[0][1], 0, 1, colors='r', linestyles='dashed')

    ax[1].plot(range(grid_num+1), np.insert(np.cumsum(mass_1d), 0, 0))
    plt.sca(ax[1])
    plt.xticks(tick_positions_m1D, tick_labels_m, size=tick_size)
    plt.yticks(size=tick_size)
    plt.xlabel(r'Companion mass ($M_{Jup}$)', size = label_size)
    plt.title('Mass CDF', size=title_size)

    ax[1].hlines(0, 0, grid_num-1, colors='k', linestyles='solid')
    ax[1].vlines(twosig_inds[1][0], 0, 1, colors='r', linestyles='dashed')
    ax[1].vlines(twosig_inds[1][1], 0, 1, colors='r', linestyles='dashed')
    
    save_dir = os.path.join(outdir, 'results/{}/'.format(star_name)) # Each star gets its own folder
    os.makedirs(save_dir, exist_ok = True)
    # if not os.path.isdir(save_dir):
    #     os.makedirs(save_dir)
    
    fig.tight_layout()
    fig.savefig(save_dir + star_name + '_1d.png')
    
    return


def tick_function_a(sep, d_star):
    """
    Converts separation in AU into
    separation in arcsec
    
    Arguments:
        sep (float, AU): Projected physical separation of
                         a companion from its host
        d_star (float, AU): Distance of host star from Earth
    
    Returns:
        ang_sep (float, arcsec): Projected angular host-
                                 companion separation
    """
    pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
    
    d_star_pc = np.array(d_star)/np.array(pc_in_au)
    
    ang_sep = sep/d_star_pc
    

    return ang_sep
        
        
        
        
        