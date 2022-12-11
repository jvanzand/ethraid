import os
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.patches as ptch

from trends import helper_functions_general as hlp


def scatter_companion(scatter_plot, grid_num_2d, a_lim_plot, m_lim_plot):
    
    sma, mass = scatter_plot
    
    sep_ind = hlp.value2index(sma, (0, grid_num_2d-1), a_lim_plot)
    mp_ind  = hlp.value2index(mass, (0, grid_num_2d-1), m_lim_plot)
    
    return sep_ind, mp_ind
        

def period_lines(m_star, a_lim, m_lim, a_lim_plot, m_lim_plot, grid_num_2d, n, how='tot'):
    
    
    a_min, m_min = a_lim[0], m_lim[0]
    a_max, m_max = a_lim[1], m_lim[1]
    
    ######## Adding lines of constant period ##########
    hip_times  = [Time(1989.85, format='decimalyear').jd, Time(1993.21, format='decimalyear').jd]       
    #https://www.cosmos.esa.int/web/hipparcos/catalogue-summary

    gaia_times = [Time('2014-07-25', format='isot').jd, Time('2017-05-28', format='isot').jd] 
    #https://www.cosmos.esa.int/web/gaia/earlydr3

    # Time between the midpoints of the two missions
    baseline_days = ((gaia_times[1] + gaia_times[0])/2 - (hip_times[1] + hip_times[0])/2)
    gaia_baseline_days = gaia_times[1] - gaia_times[0]

    # Log-spaced masses in Jupiter masses
    const_per_m_list = np.logspace(np.log10(m_min), np.log10(m_lim[1]), 50)
    const_per_m_inds = hlp.value2index(const_per_m_list, (0, grid_num_2d-1), m_lim_plot)
    
    if how == 'tot':
        baseline = baseline_days
        fmt = '--k'
    else:
        baseline = gaia_baseline_days
        fmt = '--r'

    # Lines of constant period for p = baseline_days/n

    const_per_a_list = hlp.period_lines(const_per_m_list, baseline/(n+1), m_star)
    const_per_a_inds = hlp.value2index(const_per_a_list, (0, grid_num_2d-1), a_lim_plot)

        
    values_in_bounds = np.where((a_lim[0] < const_per_a_list)&(const_per_a_list < a_lim[1]))
    
    return const_per_a_inds[values_in_bounds], const_per_m_inds[values_in_bounds], fmt

        
def marginalized_1d(star_name, post_tot, grid_num, twosig_inds, a_lim, m_lim, tick_labels_a, tick_labels_m):
    
    title_size = 30
    label_size = 25
    tick_num = 6
    tick_size = 25
    
    fig, ax = plt.subplots(1,2, figsize=(12,8))
    sma_1d = post_tot.sum(axis=0)
    mass_1d = post_tot.sum(axis=1)
    
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
    
    save_dir_1D = 'results/1D_posts/' # 2D images of all stars in one folder, 1D images in another
    if not os.path.isdir(save_dir_1D):
        os.makedirs(save_dir_1D)
    
    fig.tight_layout()
    fig.savefig(save_dir_1D + star_name + '_1d.png')
    
    return


def tick_function_a(sep, d_star):
    """
    Converts separation in AU into
    separation in arcsec
    """
    pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
    
    d_star_pc = np.array(d_star)/np.array(pc_in_au)
    

    return sep/d_star_pc  
        
        
        
        
        