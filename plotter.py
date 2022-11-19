import os
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.patches as ptch

import helper_functions_general as hlp
import helper_functions_plotting as hlp_plot

def joint_plot(star_name, m_star, d_star, vmag, post_tot, post_rv, post_astro, post_imag, grid_num, a_lim, m_lim,
               scatter_plot=None, period_lines=False, marginalized=True):
    
    tick_num = 6
    tick_size = 40
    
    a_min, m_min = a_lim[0], m_lim[0]
    a_max, m_max = a_lim[1], m_lim[1]
    
    fig, ax = plt.subplots(figsize=(12,12), dpi = 300)
    
    ######## Padding arrays #########
    
    grid_pad = int(np.round(grid_num/15)) # grid_pad is the number of index blocks by which the grid is padded
    
    frac_exp = grid_pad/grid_num # This is the fraction by which the grid is extended to include the ruled out regions. Since it's log scale, this is an exponent.
    
    
    a_min_plot = a_min/(a_max/a_min)**(frac_exp)
    m_min_plot = m_min/(m_max/m_min)**(frac_exp)
    
    
    ######################
    # import helper_functions_imaging as hlp_imag
    # post_imag2 = hlp_imag.imag_array(d_star, vmag, 3.77, 'data/EDG_clean_curves/vortex_Lband.csv', a_lim, m_lim, 100)
    # post_imag2_pad = np.pad(post_imag2, [(grid_pad, 0), (grid_pad, 0)])
    # t_contours_imag2 = hlp.contour_levels(post_imag2, [1,2])
    # post_imag_cont2 = ax.contourf(post_imag2_pad, t_contours_imag2,
    #                            cmap='gray', extend='max', alpha=0.7, zorder=1)
    ######################
    
    post_imag_pad = np.pad(post_imag, [(grid_pad, 0), (grid_pad, 0)])
    post_rv_pad = np.pad(post_rv, [(grid_pad, 0), (grid_pad, 0)])
    post_astro_pad = np.pad(post_astro, [(grid_pad, 0), (grid_pad, 0)])
    post_tot_pad = np.pad(post_tot, [(grid_pad, 0), (grid_pad, 0)])
    
    try:
        t_contours_astro = hlp.contour_levels(post_astro, [1,2])
        post_astro_cont = ax.contourf(post_astro_pad, t_contours_astro,
                         cmap='Blues', extend='max', alpha=0.5, zorder=10)
    
    except:
        print('Error encountered in astrometry plot. Moving on.')
        pass
    
    t_contours_imag = hlp.contour_levels(post_imag, [1,2])
    t_contours_rv = hlp.contour_levels(post_rv, [1,2])
    t_contours_tot = hlp.contour_levels(post_tot, [1,2])
    
    # Un-filled contour for imaging to just show a line
    post_imag_cont = ax.contour(post_imag_pad, t_contours_imag,
                               cmap='gray', extend='max', alpha=0.4, zorder=0)
    post_rv_cont = ax.contourf(post_rv_pad, t_contours_rv,
                               cmap='Greens', extend='max', alpha=0.5, zorder=20)
    post_tot_cont = ax.contourf(post_tot_pad, t_contours_tot,
                               cmap='Reds', extend='max', alpha=0.75, zorder=30)
    
    # grid_num_2d is the side length of the 2D plotting array
    grid_num_2d = grid_num+grid_pad
    
    # We want the rectangles to be grid_num_2d long, and grid_pad wide
    mass_rect = ptch.Rectangle((0, 0), grid_num_2d-1, grid_pad,
                                       color='gray', alpha=1.0)
    a_rect = ptch.Rectangle((0, 0), grid_pad, grid_num_2d-1,
                                       color='gray', alpha=1.0)

    ax.add_patch(mass_rect)
    ax.add_patch(a_rect)
    
    ############### In-plot Labels #####################
    label_size = 50
    region_label_size = 50
    restricted_region_label_size = 40

    plt.text((5/16)*grid_num_2d, (1/8)*(grid_pad/2), 'Ruled out by RVs', 
              size=restricted_region_label_size)
    plt.text((1/4)*(grid_pad/2), (1/16)*grid_num_2d, 'Ruled out by minimum period', 
              size=restricted_region_label_size, rotation=90)


    ax.set_xlabel('Semi-major axis (au)', size=label_size)
    ax.set_ylabel(r'$M_p$ ($M_{Jup}$)', size=label_size)
    # ax.set_ylabel(r'$M_p$ ($M_{\odot}$)', size=label_size)
    ###################################################
    ############ Axis ticks and labels ################
    
    # List of round numbers to use as labels for both a and m
    # tick_labels = np.array([0.11, 0.33, 1, 3, 10, 30, 100, 300, 900])
    min_exp = -4
    max_exp = 13
    n = max_exp-min_exp+1
    # tick_labels = np.array([0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])
    tick_labels = np.logspace(min_exp, max_exp, n, base=4)

    # Chop out any labels outside the a or m bounds
    raw_labels_a = tick_labels[(a_lim[0] < tick_labels) & (tick_labels < a_lim[1])][:tick_num]
    raw_labels_m = tick_labels[(m_lim[0] < tick_labels) & (tick_labels < m_lim[1])][:tick_num]


    # Make sure the whole numbers are integers for clean display, but the small floats are rounded to 2 decimals
    tick_labels_a = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), raw_labels_a))
    tick_labels_m = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), raw_labels_m))
    
    
    # Convert the labels to index positions. Note that the positions need not be integers, even though they correspond to "indices"
    a_lim_plot = (a_min_plot, a_max)
    m_lim_plot = (m_min_plot, m_max)
    
    tick_positions_a = hlp.value2index(tick_labels_a, (0, grid_num_2d-1), a_lim_plot)
    tick_positions_m = hlp.value2index(tick_labels_m, (0, grid_num_2d-1), m_lim_plot)
    
    
    # #############################
    # M_sun = 1.988409870698051e+33
    # M_jup = 1.8981245973360504e+30
    #
    # raw_labels_msun = np.array([0.01, 0.05, 0.25, 1.0])
    # raw_labels_mjup = raw_labels_msun*M_sun/M_jup
    #
    # tick_labels_m = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), raw_labels_msun))
    # tick_positions_m = hlp.value2index(raw_labels_mjup, (0, grid_num_2d-1), m_lim_plot)

    # #############################
    
    plt.xticks(tick_positions_a, tick_labels_a, size=tick_size)
    plt.yticks(tick_positions_m, tick_labels_m, size=tick_size)
    
    
    ######## Done with x and y axes. Now to add the top x axis, which is separation in arcseconds
    raw_labels_sep = hlp_plot.tick_function_a(tick_labels_a, d_star)
    tick_labels_sep = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), raw_labels_sep))

    ax2 = ax.twiny()
    plt.sca(ax2)
    plt.xlim(0, grid_num+grid_pad-1)
    plt.xticks(tick_positions_a, tick_labels_sep, size=tick_size*0.75)
    plt.xlabel('Angular separation (arcsec)', size=label_size*0.75)
    
    
    
    
    if scatter_plot is not None:
        
        sep_ind, mp_ind  = hlp_plot.scatter_companion(scatter_plot, grid_num_2d, a_lim_plot, m_lim_plot)

        plt.scatter(sep_ind, mp_ind, marker='*', c='yellow', edgecolors='black', s=2000, zorder=4)
    
    if period_lines:
        
        for n in range(5):
            const_per_a_inds, const_per_m_inds, fmt =\
                                    hlp_plot.period_lines(m_star, a_lim, m_lim, 
                                                          a_lim_plot, m_lim_plot, 
                                                          grid_num_2d, n, how='tot')
            plt.plot(const_per_a_inds, const_per_m_inds, fmt, alpha=0.5)
            
            
            const_per_a_inds, const_per_m_inds, fmt =\
                                    hlp_plot.period_lines(m_star, a_lim, m_lim, 
                                                          a_lim_plot, m_lim_plot, 
                                                          grid_num_2d, n, how='gaia')                             
            plt.plot(const_per_a_inds, const_per_m_inds, fmt, alpha=0.5)
            


    fig.tight_layout()
    # save_dir_2D = 'results/'+star_name+'/' # Each star gets its own folder
    save_dir_2D = 'results/2D_posts/' # 2D images of all stars in one folder, 1D images in another
    # Try to make directory. If it exists, just continue. Parallel code was bugging out here, so exist_ok is great.
    os.makedirs(save_dir_2D, exist_ok = True)
    fig.savefig(save_dir_2D + star_name + '.png')
    plt.close()
    
    # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # twosig_inds contains the indices where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)

    if marginalized:
        
        hlp_plot.marginalized_1d(star_name, post_tot, grid_num, twosig_inds, 
                                 a_lim, m_lim, tick_labels_a, tick_labels_m)

        
    # Print out the 2-sigma boundaries (bounds) for the joint posterior
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained in the part of the posterior inside of which a horizontal line at height two_sig_levels[i] falls.
    print('a_lim = ', bounds[0], ' AU')
    print('m_lim = ', bounds[1], ' M_J')
    
    return
