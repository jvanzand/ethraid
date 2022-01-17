import os
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.patches as ptch

import helper_functions_general as hlp


def joint_plot(star_name, m_star, post_tot, post_rv, post_astro, grid_num, a_lim, m_lim,
               period_lines = False, marginalized=True):
    
    tick_num = 6
    tick_size = 30
    
    a_min, m_min = a_lim[0], m_lim[0]
    a_max, m_max = a_lim[1], m_lim[1]
    
    fig, ax = plt.subplots(figsize=(12,12), dpi = 300)
    
    ######## Padding arrays #########
    
    grid_pad = int(np.round(grid_num/15)) # grid_pad is the number of index blocks by which the grid is padded
    
    frac_exp = grid_pad/grid_num # This is the fraction by which the grid is extended. Since it's log scale, this is an exponent
    
    
    a_min_plot = a_min/(a_max/a_min)**(frac_exp)
    m_min_plot = m_min/(m_max/m_min)**(frac_exp)
    
    
    post_rv_pad = np.pad(post_rv, [(grid_pad, 0), (grid_pad, 0)])
    post_astro_pad = np.pad(post_astro, [(grid_pad, 0), (grid_pad, 0)])
    post_tot_pad = np.pad(post_tot, [(grid_pad, 0), (grid_pad, 0)])
    
    
    try:
        t_contours_astro = hlp.contour_levels(post_astro, [1,2]) ## !! MAybe change this to post_astro_pad
        post_astro_cont = ax.contourf(post_astro_pad, t_contours_astro, cmap='Blues', extend='max', alpha=0.5)
    
    except:
        pass
    
    t_contours_rv = hlp.contour_levels(post_rv, [1,2])
    t_contours_tot = hlp.contour_levels(post_tot, [1,2])

    post_rv_cont = ax.contourf(post_rv_pad, t_contours_rv, 
                               cmap='Greens', extend='max', alpha=0.5)
    post_tot_cont = ax.contourf(post_tot_pad, t_contours_tot,
                       cmap='Reds', extend='max', alpha=0.75)

    
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
    restricted_region_label_size = 35

    plt.text((5/16)*grid_num_2d, (1/4)*(grid_pad/2), 'Ruled out by RVs', 
              size=restricted_region_label_size)
    plt.text((1/4)*(grid_pad/2), (1/8)*grid_num_2d, 'Ruled out by minimum period', 
              size=restricted_region_label_size, rotation=90)


    ax.set_xlabel('Semi-major Axis (au)', size=label_size)
    ax.set_ylabel(r'$M_p$ ($M_{Jup}$)', size=label_size)
    ###################################################
    ############ Axis ticks and labels ################
    
    # List of round numbers to use as labels for both a and m
    #tick_labels = np.array([0.11, 0.33, 1, 3, 10, 30, 100, 300, 900])
    tick_labels = np.array([0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024])

    # Chop out any labels outside the a or m bounds
    tick_labels_a = tick_labels[(a_lim[0] < tick_labels) & (tick_labels < a_lim[1])]
    tick_labels_m = tick_labels[(m_lim[0] < tick_labels) & (tick_labels < m_lim[1])]

    # Make sure the whole numbers are integers for clean display, but the small floats are rounded to 2 decimals
    tick_labels_a = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), tick_labels_a))
    tick_labels_m = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), tick_labels_m))
    
    
    # Convert the labels to index positions. Note that the positions need not be integers, even though they correspond to "indices"
    a_lim_plot = (a_min_plot, a_max)
    m_lim_plot = (m_min_plot, m_max)
    
    tick_positions_a = hlp.value2index(tick_labels_a, (0, grid_num_2d-1), a_lim_plot)
    tick_positions_m = hlp.value2index(tick_labels_m, (0, grid_num_2d-1), m_lim_plot)
    
    plt.xticks(tick_positions_a, [str(i) for i in tick_labels_a], size=tick_size)
    plt.yticks(tick_positions_m, [str(i) for i in tick_labels_m], size=tick_size)
    
    ### Experimental: adding top x-axis to show separations
    # def au2sep(au):
    #     """
    #     Given system distance in pc, converts separation in au into separation in arcsec
    #     """
    #     asec = au/10
    #     return asec
    #
    # def sep2au(asec):
    #     """
    #     Given system distance in pc, converts separation in arcsec into separation in au
    #     """
    #     au = asec*10
    #     return au
    #
    # ax.secondary_xaxis('top', functions=(au2sep, sep2au))
    
    
    if period_lines:
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

        # Lines of constant period for p = baseline_days/n
        for f in range(5):

            const_per_a_list = hlp.period_lines(const_per_m_list, baseline_days/(f+1), m_star)
            const_per_a_inds = hlp.value2index(const_per_a_list, (0, grid_num_2d-1), a_lim_plot)

            
            values_in_bounds = np.where((a_lim[0] < const_per_a_list)&(const_per_a_list < a_lim[1]))
            

            plt.plot(const_per_a_inds[values_in_bounds], const_per_m_inds[values_in_bounds], '--k', alpha=0.5)
            #plt.plot(const_per_a_inds, const_per_m_inds, '--k', alpha=0.5)

        # Lines of constant period for p = gaia_baseline_days/n
        for f in range(5):

            const_per_a_list = hlp.period_lines(const_per_m_list, gaia_baseline_days/(f+1), m_star)
            const_per_a_inds = hlp.value2index(const_per_a_list, (0, grid_num_2d-1), a_lim_plot)

            values_in_bounds = np.where(const_per_a_list >= a_min)

            plt.plot(const_per_a_inds[values_in_bounds], const_per_m_inds[values_in_bounds], '--r', alpha=0.5)


    fig.tight_layout()
    # save_dir_2D = 'results/'+star_name+'/' # Each star gets its own folder
    save_dir_2D = 'results/2D_posts/' # 2D images of all stars in one folder, 1D images in another
    if not os.path.isdir(save_dir_2D):
        os.makedirs(save_dir_2D)
    fig.savefig(save_dir_2D + star_name + '.png')
    plt.close()
    
    # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # twosig_inds contains the indices where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)

    if marginalized:
        
        title_size = 30
        label_size = 25
        tick_num = 6
        tick_size = 25
        
        fig, ax = plt.subplots(1,2, figsize=(12,8))
        sma_1d = post_tot.sum(axis=0)
        mass_1d = post_tot.sum(axis=1)

        ax[0].plot(range(grid_num+1), np.insert(np.cumsum(sma_1d), 0, 0))
        plt.sca(ax[0])
        plt.xticks(tick_positions_a, tick_labels_a, size=tick_size)
        plt.yticks(size=tick_size)
        plt.title('Semi-major axis CDF', size=title_size)
        plt.xlabel('Companion semi-major axis (AU)', size = label_size)

        ax[0].hlines(0, 0, grid_num-1, colors='k', linestyles='solid')
        ax[0].vlines(twosig_inds[0][0], 0, 1, colors='r', linestyles='dashed')
        ax[0].vlines(twosig_inds[0][1], 0, 1, colors='r', linestyles='dashed')

        ax[1].plot(range(grid_num+1), np.insert(np.cumsum(mass_1d), 0, 0))
        plt.sca(ax[1])
        plt.xticks(tick_positions_m, tick_labels_m, size=tick_size)
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
        
    # Print out the 2-sigma boundaries (bounds) for the joint posterior
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained in the part of the posterior inside of which a horizontal line at height two_sig_levels[i] falls.
    print('a_lim = ', bounds[0], ' AU')
    print('m_lim = ', bounds[1], ' M_J')
    
    return
