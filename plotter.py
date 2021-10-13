import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.patches as ptch

from trends import helper_functions_wrapper as hlpw


def joint_plot(m_star, post_tot, post_rv, post_astro, grid_num, a_lim, m_lim, 
                min_vals, save_name='companion', period_lines = False, marginalized=True):
    
    tick_num = 6
    tick_size = 30
    
    min_a, min_m = min_vals
    
    fig, ax = plt.subplots(figsize=(12,12), dpi = 300)
    
    try:
        t_contours_astro = hlpw.contour_levels(post_astro, [1,2])
        post_astro_cont = ax.contourf(post_astro, t_contours_astro, cmap='Blues', extend='max', alpha=0.5)
    
    except:
        pass
    
    t_contours_rv = hlpw.contour_levels(post_rv, [1,2])

    post_rv_cont = ax.contourf(post_rv, t_contours_rv,
                               cmap='Greens', extend='max', alpha=0.5)
                               
    t_contours_tot = hlpw.contour_levels(post_tot, [1,2])
    post_tot_cont = ax.contourf(post_tot, t_contours_tot,
                       cmap='Reds', extend='max', alpha=0.75)

    a_list = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
    m_list = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)

    min_index_m = hlpw.value2index(min_m, (0, grid_num-1), m_lim)
    min_index_a = hlpw.value2index(min_a, (0, grid_num-1), a_lim)


    mass_rect = ptch.Rectangle((0, 0), grid_num-1, min_index_m, 
                                       color='gray', alpha=1.0)
    a_rect = ptch.Rectangle((0, 0), min_index_a, grid_num-1, 
                                       color='gray', alpha=1.0)

    ax.add_patch(mass_rect)
    ax.add_patch(a_rect)
    
    ############### In-plot Labels #####################
    label_size = 50
    region_label_size = 50
    restricted_region_label_size = 35

    # plt.text((16/32)*grid_num, (7/8)*grid_num, 'RV',
    #           size=region_label_size, rotation=50)
    # plt.text((7/16)*grid_num, (1/4)*grid_num, 'Astrometry',
    #           size=region_label_size)

    plt.text((1/6)*grid_num, (1/3)*(min_index_m-1), 'Masses disallowed by RVs', 
              size=restricted_region_label_size)
    plt.text((1/3)*(min_index_a-1), (1/8)*grid_num, 'Ruled out by minimum period', 
              size=restricted_region_label_size, rotation=90)


    ax.set_xlabel('Semi-major Axis (au)', size=label_size)
    ax.set_ylabel(r'$M_p$ ($M_{Jup}$)', size=label_size)
    ###################################################
    ############ Axis ticks and labels ################
    
    # List of round numbers to use as labels for both a and m
    #tick_labels = np.array([0.11, 0.33, 1, 3, 10, 30, 100, 300, 900, 2700])
    tick_labels = np.array([0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048])

    # Chop out any labels outside the a or m bounds
    tick_labels_a = tick_labels[(a_lim[0] < tick_labels) & (tick_labels < a_lim[1])]
    tick_labels_m = tick_labels[(m_lim[0] < tick_labels) & (tick_labels < m_lim[1])]

    # Make sure the whole numbers are integers for clean display, but the small floats are rounded to 2 decimals
    tick_labels_a = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), tick_labels_a))
    tick_labels_m = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), tick_labels_m))

    # Convert the labels to index positions. Note that the positions need not be integers, even though they correspond to "indices"
    tick_positions_a = hlpw.value2index(tick_labels_a, (0, grid_num-1), a_lim)
    tick_positions_m = hlpw.value2index(tick_labels_m, (0, grid_num-1), m_lim)
    
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
        const_per_m_list = np.logspace(np.log10(min_m), np.log10(m_lim[1]), 50)
        const_per_m_inds = hlpw.value2index(const_per_m_list, (0, grid_num-1), m_lim)

        # Lines of constant period for p = baseline_days/n
        for f in range(5):

            const_per_a_list = hlpw.period_lines(const_per_m_list, baseline_days/(f+1), m_star)
            const_per_a_inds = hlpw.value2index(const_per_a_list, (0, grid_num-1), a_lim)

            
            values_in_bounds = np.where((a_lim[0] < const_per_a_list)&(const_per_a_list < a_lim[1]))
            

            plt.plot(const_per_a_inds[values_in_bounds], const_per_m_inds[values_in_bounds], '--k', alpha=0.5)
            #plt.plot(const_per_a_inds, const_per_m_inds, '--k', alpha=0.5)

        # Lines of constant period for p = gaia_baseline_days/n
        for f in range(5):

            const_per_a_list = hlpw.period_lines(const_per_m_list, gaia_baseline_days/(f+1), m_star)
            const_per_a_inds = hlpw.value2index(const_per_a_list, (0, grid_num-1), a_lim)

            values_in_bounds = np.where(const_per_a_list >= min_a)

            plt.plot(const_per_a_inds[values_in_bounds], const_per_m_inds[values_in_bounds], '--r', alpha=0.5)
    
    ###########
    synth_a = 3.5
    synth_m = 15

    inds = hlpw.value2index(synth_a, (0, grid_num-1), a_lim),\
           hlpw.value2index(synth_m, (0, grid_num-1), m_lim)

    ax.scatter(inds[0], inds[1], s=600, c='yellow', marker='*')
    ###########

    fig.tight_layout()
    fig.savefig('plots/' + save_name + '.png')
    plt.close()
    
    # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained in the part of the posterior inside of which a horizontal line at height two_sig_levels[i] falls.
    # twosig_inds contains the indices where the above horizontal line crosses the posterior. In case it crosses more than twice, it contains the first and last instances.
    bounds, twosig_levels, twosig_inds = hlpw.bounds_1D(post_tot, [m_lim, a_lim], interp_num = 1e4)
    
    
    if marginalized:
        
        title_size = 30
        label_size = 25
        tick_num = 6
        tick_size = 25
        
        fig, ax = plt.subplots(1,2, figsize=(12,8))
        sma_1d = post_tot.sum(axis=0)
        mass_1d = post_tot.sum(axis=1)

        ax[0].plot(range(grid_num), sma_1d)
        plt.sca(ax[0])
        plt.xticks(tick_positions_a, tick_labels_a, size=tick_size)
        plt.yticks(size=tick_size)
        plt.title('Semi-major axis posterior', size=title_size)
        plt.xlabel('Companion semi-major axis (AU)', size = label_size)
        ax[0].hlines(twosig_levels[0], 0, grid_num-1, colors='k', linestyles='solid')
        ax[0].vlines(twosig_inds[0][0], 0, 3*twosig_levels[0], colors='r', linestyles='dashed')
        ax[0].vlines(twosig_inds[0][1], 0, 3*twosig_levels[0], colors='r', linestyles='dashed')

        ax[1].plot(range(grid_num), mass_1d)
        plt.sca(ax[1])
        plt.xticks(tick_positions_m, tick_labels_m, size=tick_size)
        plt.yticks(size=tick_size)
        plt.xlabel(r'Companion mass ($M_{Jup}$)', size = label_size)
        plt.title('Mass posterior', size=title_size)
        ax[1].hlines(twosig_levels[1], 0, grid_num-1, colors='k', linestyles='solid')
        ax[1].vlines(twosig_inds[1][0], 0, 3*twosig_levels[1], colors='r', linestyles='dashed')
        ax[1].vlines(twosig_inds[1][1], 0, 3*twosig_levels[1], colors='r', linestyles='dashed')
        
        fig.tight_layout()
        fig.savefig('plots/' + save_name + '_1d.png')
        
    # Print out the 2-sigma boundaries (bounds) for the joint posterior
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained in the part of the posterior inside of which a horizontal line at height two_sig_levels[i] falls.
    print('a_lim = ', bounds[0], ' AU')
    print('m_lim = ', bounds[1], ' M_J')
    
    return
