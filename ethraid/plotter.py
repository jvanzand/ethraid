import os
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import warnings

from ethraid import helper_functions_plotting as hlp_plot
from ethraid.compiled import helper_functions_general as hlp

def joint_plot(star_name, m_star, d_star, 
               run_rv, run_astro, run_imag, 
               post_tot, post_rv, post_astro, post_imag, 
               grid_num, a_lim, m_lim,
               scatter_plot=None, period_lines=False, 
               outdir='', verbose=False):
    
    """
    Plots the 2D joint mass-semi-major axis posteriors calculated using 
    provided RV, astrometry, and imaging data.

    Arguments:
        star_name (str): Name of star (does not need to be official)
        m_star (float, M_jup): Mass of host star
        d_star (float, AU): Distance from Earth to host star
               
        run_rv (bool): Was RV data used in calculation?
        run_astro (bool): Was astrometry data used in calculation?
        run_imag (bool): Was imaging data used in calculation?
               
        post_tot (array of floats): Total posterior array
        post_rv (array of floats): Model probabilities given RV data only
        post_astro (array of floats): Model probabilities given astrometry data only
        post_imag (array of floats): Model probabilities given imaging data only
               
        grid_num (int): Shape of square posterior arrays
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max)
        m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
        run_rv, run_astro, run_imag (bool): True/False values indicating whether
                                            each data type was considered. Omitted
                                            data types are not plotted.
        scatter_plot (list of floats): Optional (semi-major axis, mass) pair
                                        specifying the location of a known
                                        companion to plot. Sma in AU, mass in
                                        M_jup.
         period_lines (bool): Optionally plot lines of constant period
                              at periods equal to harmonics of the Gaia and 
                              HG baselines
         out_dir (str): Path to save generated plot

    Returns:
         None (plots 2D joint posterior)
    """
    
    a_min, m_min = a_lim[0], m_lim[0]
    a_max, m_max = a_lim[1], m_lim[1]
    
    fig, ax = plt.subplots(figsize=(12,12), dpi = 300)
    
    ######## Padding arrays #########
    ## Companions of low mass and close separation are ruled out based on the RV trend alone. To illustrate this, expand the grid slightly on the left and bottom.
    
    grid_pad = int(np.round(grid_num/15)) # grid_pad is the number of index blocks by which the grid is padded
    
    frac_exp = grid_pad/grid_num # This is the fraction by which the grid is extended to include the ruled out regions. Since it's log scale, this is an exponent.
      
    a_min_plot = a_min/(a_max/a_min)**(frac_exp)
    m_min_plot = m_min/(m_max/m_min)**(frac_exp)
                 
         
    # Loop through all data types (rv, astro, and imag) simultaneously define and format variables, as well as the values those variables are being assigned.
    ################################
    data_types = ['rv', 'astro', 'imag']
    colors = ['Greens', 'Blues', 'gray']
    alphas = [0.5, 0.5, 0.4]
    zorders = [20, 10, 0]
    
    
    for i in range(3):
        dt = data_types[i]
        c = colors[i]
        alpha = alphas[i]
        z = zorders[i]
        # First, check if run_rv, run_astro, and run_imag are True. If any is not, don't plot that contour
        if eval("run_{}".format(dt)):
            exec("post_{0}_pad = np.pad(post_{0}, [(grid_pad, 0), (grid_pad, 0)])".format(dt))
            
            # For imaging only, plot 2sig contour and use contour instead of contourf to get a line instead of a filled region
            if dt == 'imag':
                exec("t_contours_{0} = hlp.contour_levels(post_{0}, [1,2])".format(dt))
                
                # In the approximate case, the imaging posterior has 2 regions by design: uniformly 0 and uniformly some nonzero value. Suppress Matplotlib's warning that contour levels are undefined in this case.
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message='No contour levels were found within the data range.')
                    exec("post_{0}_cont = ax.contour(post_{0}_pad, t_contours_{0},\
                                 cmap='{1}', extend='max', alpha={2}, zorder={3})".format(dt, c, alpha, z))
                                 
            else:
                exec("t_contours_{0} = hlp.contour_levels(post_{0}, [1,2])".format(dt))

                exec("post_{0}_cont = ax.contourf(post_{0}_pad, t_contours_{0},\
                             cmap='{1}', extend='max', alpha={2}, zorder={3})".format(dt, c, alpha, z))
            
        
    post_tot_pad = np.pad(post_tot, [(grid_pad, 0), (grid_pad, 0)])
    t_contours_tot = hlp.contour_levels(post_tot, [1,2])
    post_tot_cont = ax.contourf(post_tot_pad, t_contours_tot,
                                cmap='Reds', extend='max', alpha=0.75, zorder=30)   
    ################################
    
    
    # ########
    # ### Optionally add Vortex coronagraph contrast curve to 2D plot
    # from ethraid import helper_functions_imaging as hlp_imag
    # contrast_str = 'ethraid/data/clean_curves/vortex_Lband.csv'
    # post_vortex = hlp_imag.imag_array(d_star, vmag, 3.77, contrast_str, a_lim, m_lim, grid_num)
    #
    #
    # post_vortex_pad = np.pad(post_vortex, [(grid_pad, 0), (grid_pad, 0)])
    # t_contours_vortex = hlp.contour_levels(post_vortex, [1,2])
    #
    # post_vortex_cont = ax.contour(post_vortex_pad, t_contours_vortex,
    #                            cmap='autumn', extend='both', alpha=1.0, zorder=35)
    # #######
    
    
    # grid_num_ext is the side length of the 2D plotting array
    grid_num_ext = grid_num+grid_pad
    
    # We want the rectangles to be grid_num_ext long, and grid_pad wide
    mass_rect = ptch.Rectangle((0, 0), grid_num_ext-1, grid_pad,
                                       color='gray', alpha=1.0, zorder=100)
    a_rect = ptch.Rectangle((0, 0), grid_pad, grid_num_ext-1,
                                       color='gray', alpha=1.0, zorder=100)

    ax.add_patch(mass_rect)
    ax.add_patch(a_rect)
    
    ############### In-plot Labels #####################
    label_size = 50
    region_label_size = 50
    restricted_region_label_size = 40

    # plt.text((5/16)*grid_num_ext, (1/8)*(grid_pad/2), 'Ruled out by RVs',
    #           size=restricted_region_label_size, zorder=101)
    # plt.text((1/4)*(grid_pad/2), (1/16)*grid_num_ext, 'Ruled out by minimum period',
    #           size=restricted_region_label_size, rotation=90, zorder=101)


    ax.set_xlabel('Semi-major axis (au)', size=label_size)
    ax.set_ylabel(r'$M_p$ ($M_{Jup}$)', size=label_size)

    ###################################################
    ############ Axis ticks and labels ################
    tick_num = 6
    tick_size = 40
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
    
    tick_positions_a = hlp.value2index(tick_labels_a, (0, grid_num_ext-1), a_lim_plot)
    tick_positions_m = hlp.value2index(tick_labels_m, (0, grid_num_ext-1), m_lim_plot)

    
    plt.xticks(tick_positions_a, tick_labels_a, size=tick_size)
    plt.yticks(tick_positions_m, tick_labels_m, size=tick_size)
    
    ## Removed angular separation axis. Orbits are not face-on/circular in general, so 'a' does not correspond uniquely to angular separation.
    ########
    ######## Done with x and y axes. Now to add the top x axis, which is separation in arcseconds ########
    # raw_labels_sep = hlp_plot.tick_function_a(tick_labels_a, d_star)
    # tick_labels_sep = list(map(lambda x: int(x) if x%1 == 0 else np.around(x, decimals=2), raw_labels_sep))
    #
    # ax2 = ax.twiny()
    # plt.sca(ax2)
    # plt.xlim(0, grid_num+grid_pad-1)
    # plt.xticks(tick_positions_a, tick_labels_sep, size=tick_size*0.75)
    # plt.xlabel('Angular separation (arcsec)', size=label_size*0.75)
    
    ## Add scatter point to indicate known/expected companion location ##
    if scatter_plot is not None:

        sep_ind, mp_ind  = hlp_plot.scatter_companion(scatter_plot, grid_num_ext, a_lim_plot, m_lim_plot)

        plt.scatter(sep_ind, mp_ind, marker='*', c='yellow', edgecolors='black', s=2000, zorder=40)
    
    ## Plot lines of constant period at harmonics of mission baseline (baseline/1, baseline/2, etc.)
    if period_lines:
        ## Plot harmonics of total baseline
        const_per_a_inds_list, const_per_m_inds_list, fmt =\
                                hlp_plot.period_lines(m_star, a_lim, m_lim, 
                                                      a_lim_plot, m_lim_plot, 
                                                      grid_num_ext, 3, how='tot')
        num_lines = len(const_per_a_inds_list)
        for i in range(num_lines):
            plt.plot(const_per_a_inds_list[i], const_per_m_inds_list[i], fmt, alpha=0.5)
            
        
        ## Plot harmonics of Gaia baseline
        const_per_a_inds_list, const_per_m_inds_list, fmt =\
                                hlp_plot.period_lines(m_star, a_lim, m_lim,
                                                      a_lim_plot, m_lim_plot,
                                                      grid_num_ext, 3, how='gaia')
        for i in range(num_lines):
            plt.plot(const_per_a_inds_list[i], const_per_m_inds_list[i], fmt, alpha=0.5)
            


    fig.tight_layout()
    save_dir = os.path.join(outdir, 'results/{}/'.format(star_name)) # Each star gets its own folder
    # Try to make directory. If it exists, just continue. Parallel code was bugging out here, so exist_ok is great.
    os.makedirs(save_dir, exist_ok = True)
    fig.savefig(save_dir + star_name + '_2d.png')
    plt.close()
    
    # # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # # twosig_inds contains the indices where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    # bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)
    
    return


def plot_1d(star_name, post_tot, a_lim, m_lim, outdir=''):
    
    """
    Plots the 1D joint mass and semi-major axis posteriors calculated with 
    provided RV, astrometry, and imaging data.

    Arguments:
        star_name (str): Name of star (does not need to be official)
        post_tot (array of floats): Total posterior array
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max)
        m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
        out_dir (str): Path to save generated plot

    Returns:
         None (plots 1D posteriors)
    """
               
    ###################################################
    ############ Axis ticks and labels ################
    tick_num = 6
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
    
    ###################################################
    ###################################################
    
    # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # twosig_inds contains the indices corresponding to bounds. That is, where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)


    hlp_plot.marginalized_1d(star_name, post_tot, twosig_inds, 
                             a_lim, m_lim, tick_labels_a, tick_labels_m, 
                             outdir=outdir)
    
    return