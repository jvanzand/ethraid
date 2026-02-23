import os
import numpy as np
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.patches as ptch
import warnings

from ethraid import _ROOT
from ethraid import helper_functions_plotting as hlp_plot
from ethraid.compiled import helper_functions_general as hlp
from ethraid.compiled import helper_functions_imaging as hlp_imag

# Use plot template
plt.style.use(os.path.join(_ROOT, 'data/matplotlibrc'))

def joint_plot(star_name, m_star, d_star, 
               run_rv, run_astro, run_imag, 
               post_tot, post_rv, post_astro, post_imag, 
               a_lim, m_lim,
               scatter_plot=None, age_table=4,
               period_lines=False, outdir='', 
               verbose=False):
    
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
        scatter_plot (list): Optional list of (sep, mass) tuples to scatter plot 
                             the parameters of 1 or more companions. Sma in AU, 
                             mass in M_jup.
        age_table (int): Integer 1-5, indicating which BD cooling model to use
                         based on age of system.
                         1-->0.1 Gyr, 2-->0.5 Gyr, 3-->1 Gyr, 4-->5 Gyr, 5-->10 Gyr
                         Only needed if plotting scatter companion with angsep_mag units
        period_lines (bool): Optionally plot lines of constant period
                              at periods equal to harmonics of the Gaia and 
                              HG baselines
        out_dir (str): Path to save generated plot

    Returns:
         None (plots 2D joint posterior)
    """
    grid_num = np.shape(post_tot)[0] # grid_num is the shape of the posterior(s). This handles raw and processed posterior arrays.
    
    fig, ax = plt.subplots(figsize=(12,12), dpi = 300)
                 
         
    # Loop through all data types (rv, astro, and imag) simultaneously define and format variables, as well as the values those variables are being assigned.
    ################################
    run_list = [run_rv, run_astro, run_imag]
    post_list = [post_rv, post_astro, post_imag]
    data_types = ['rv', 'astro', 'imag']
    colors = ['Greens', 'Blues', 'gray']
    alphas = [0.5, 0.5, 0.4]
    zorders = [20, 10, 0]

    # data_types=['astro']
    # colors = ['Blues']
    
    for i in range(len(data_types)):
        run = run_list[i]
        post = post_list[i]
        dt = data_types[i]
        c = colors[i]
        alpha = alphas[i]
        z = zorders[i]
        # First, check if run_rv, run_astro, and run_imag are True. If any is not, don't plot that contour
        if run:
            
            # For imaging only, plot 2sig contour and use contour instead of contourf to get a line instead of a filled region
            if dt == 'imag':
                t_contours = hlp.contour_levels(post, [1,2])
                
                # In the approximate case, the imaging posterior has 2 regions by design: uniformly 0 and uniformly some nonzero value. Suppress Matplotlib's warning that contour levels are undefined in this case.
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message='No contour levels were found within the data range.')
                    
                    post_contour = ax.contour(post, t_contours,\
                                              cmap=c, extend='max', alpha=alpha, zorder=z)
                                 
            else:
                t_contours = hlp.contour_levels(post, [1,2])
                post_cont = ax.contourf(post, t_contours,\
                                        cmap=c, extend='max', alpha=alpha, zorder=z)
            
    ## Only plot the overlap red if plotting both RV and astro. Otherwise let the green/blue show
    if run_rv and run_astro:
        #post_tot_pad = np.pad(post_tot, [(grid_pad, 0), (grid_pad, 0)])
        t_contours_tot = hlp.contour_levels(post_tot, [1,2])
        post_tot_cont = ax.contourf(post_tot, t_contours_tot,
           cmap='Reds', extend='max', alpha=0.75, zorder=30)
    ################################
    
    ############### In-plot Labels #####################
    label_size = 50

    ax.set_xlabel('Semi-major axis (AU)', size=label_size)
    ax.set_ylabel(r'Companion Mass ($\mathrm{M_{Jup}}$)', size=label_size)

    ###################################################
    ############ Axis ticks and labels ################
    tick_num = 10
    tick_size = 40
    # List of round numbers to use as labels for both a and m
    min_exp = -3
    max_exp = 5
    n = max_exp-min_exp+1
    

    tick_vals_a_long = np.logspace(min_exp, max_exp, n, base=10)
    tick_vals_m_long = np.logspace(min_exp, max_exp, n, base=10)
    
    # Exclude any positions that fall outside the grid domain
    tick_vals_a = tick_vals_a_long[(a_lim[0]<=tick_vals_a_long)\
                                  &(tick_vals_a_long<a_lim[1])]
    tick_vals_m = tick_vals_m_long[(m_lim[0]<=tick_vals_m_long)\
                                  &(tick_vals_m_long<m_lim[1])]
    
    
    # Convert the labels to index positions. Note that the positions need not be integers, even though they correspond to "indices"
    tick_positions_a = hlp.value2index(tick_vals_a, (0, grid_num-1), a_lim)
    tick_positions_m = hlp.value2index(tick_vals_m, (0, grid_num-1), m_lim)
    
    ax.set_xticks(tick_positions_a)
    ax.set_yticks(tick_positions_m)
    ax.set_xticklabels([rf"$10^{{{int(np.log10(v))}}}$" for v in tick_vals_a], size=tick_size)
    ax.set_yticklabels([rf"$10^{{{int(np.log10(v))}}}$" for v in tick_vals_m], size=tick_size)
    
    ## Minor ticks
    minor_tick_list_of_arrs = [np.arange(2, 10) * 10.0**p for p in range(min_exp, max_exp)]
    minor_tick_vals_a_long = np.concatenate(minor_tick_list_of_arrs)
    minor_tick_vals_m_long = np.concatenate(minor_tick_list_of_arrs)
    
    # Exclude any positions that fall outside the grid domain
    minor_tick_vals_a = minor_tick_vals_a_long[(a_lim[0]<=minor_tick_vals_a_long)\
                                              &(minor_tick_vals_a_long<a_lim[1])]
    minor_tick_vals_m = minor_tick_vals_m_long[(m_lim[0]<=minor_tick_vals_m_long)\
                                              &(minor_tick_vals_m_long<m_lim[1])]
    
    minor_tick_positions_a = hlp.value2index(minor_tick_vals_a, (0, grid_num-1), a_lim)
    minor_tick_positions_m = hlp.value2index(minor_tick_vals_m, (0, grid_num-1), m_lim)
    
    ax.set_xticks(minor_tick_positions_a, minor=True)
    ax.set_yticks(minor_tick_positions_m, minor=True)
    
    
    ## Add top x-axis for angular separation ##
    ###########################################
    ## Create a top x-axis sharing the same x scale
    ax_top = ax.twiny()

    ## Match limits so ticks line up
    ax_top.set_xlim(ax.get_xlim())
    ax_top.tick_params(axis='x', which='major', pad=2)
    
    ## Initial long list of ang seps
    min_exp, max_exp = -4, 4
    n = max_exp-min_exp+1
    tick_vals_angsep_long = np.logspace(min_exp, max_exp, n, base=10)
    
    ## Max and min angular separation values, based on SMA limits
    angsep_min = hlp_plot.sma2angsep(a_lim[0], d_star)
    angsep_max = hlp_plot.sma2angsep(a_lim[1], d_star)
    angsep_lim = (angsep_min, angsep_max)
    
    ## Trim to fit in SMA domain
    tick_vals_angsep = tick_vals_angsep_long[(angsep_min<=tick_vals_angsep_long)\
                                            &(tick_vals_angsep_long<angsep_max)]
                                            
    ## Convert angsep to SMA and then to index
    tick_positions_angsep_sma = hlp_plot.angsep2sma(tick_vals_angsep, d_star)
    tick_positions_angsep = hlp.value2index(tick_positions_angsep_sma, (0, grid_num-1), a_lim)

    ax_top.set_xticks(tick_positions_angsep)
    ax_top.set_xticklabels([rf"$10^{{{int(np.log10(v))}}}$" for v in tick_vals_angsep], size=tick_size)
    
    ## Minor ticks on top axis
    minor_tick_list_of_arrs = [np.arange(2, 10) * 10.0**p for p in range(min_exp, max_exp)]
    minor_tick_vals_angsep_long = np.concatenate(minor_tick_list_of_arrs)
    minor_tick_vals_angsep = minor_tick_vals_angsep_long[(angsep_min<=minor_tick_vals_angsep_long)\
                                              &(minor_tick_vals_angsep_long<angsep_max)]
    
    minor_tick_positions_angsep = hlp.value2index(minor_tick_vals_angsep, (0, grid_num-1), angsep_lim)
    ax_top.set_xticks(minor_tick_positions_angsep, minor=True)

    ## Set label
    ax_top.set_xlabel("Projected Separation (\")", size=label_size, labelpad=15)
    #import pdb; pdb.set_trace()
    ###########################################
    
    
    ## Add scatter point to indicate known/expected companion location ##
    if scatter_plot is not None:

        for scatter_tuple in scatter_plot:
            
            ## scatter_tuple looks like (xcoord, ycoord, units)
            ## e.g. (5.1, 9.6, 'angsep_mag')
            ## For backward compatibility, default to assuming units=='sma_mass'
            try:
                units = scatter_tuple[2]
            except:
                units = 'sma_mass'
            
            
            if units=='sma_mass': # If sma and mass given directly, no need to transform
                scatter_pair = scatter_tuple[0:2]
                
            elif units=='angsep_mag': # If angsep and mag given, convert to sma and mass
                try:
                    band_name = scatter_tuple[3]
                except:
                    raise Exception("plotter.joint_plot: \n"
                                    "        If using 'angsep_mag' units, must provide \n"
                                    "        band name of mag as 4th tuple value in scatter_plot \n"
                                    "        parameter in the setup file.")
                angsep, mag = scatter_tuple[0:2]
                sma = hlp_plot.angsep2sma(angsep, d_star) # Convert ang sep to distance (ignoring projection effects)
                
                ## Now convert mag to mass
                mass = hlp_imag.mag2mass(d_star, mag, band_name, age_table)
                
                scatter_pair = (sma, mass)
                
                
            else:
                raise Exception("plotter.joint_plot: \n"
                                "        units must be either 'sma_mass' or 'angsep_mag'")
                
            sep_ind, mp_ind  = hlp_plot.scatter_companion(scatter_pair, grid_num, a_lim, m_lim)

            plt.scatter(sep_ind, mp_ind, marker='*', c='yellow', edgecolors='black', s=2000, zorder=40)
    
    ## Plot lines of constant period at harmonics of mission baseline (baseline/1, baseline/2, etc.)
    if period_lines:
        ## Plot harmonics of total baseline
        const_per_a_inds_list, const_per_m_inds_list, fmt =\
                                hlp_plot.period_lines(m_star, a_lim, m_lim, 
                                                      a_lim, m_lim, 
                                                      grid_num, 0, how='tot')
        num_lines = len(const_per_a_inds_list)
        for i in range(num_lines):
            plt.plot(const_per_a_inds_list[i], const_per_m_inds_list[i], fmt, alpha=0.5)
            
        
        ## Plot harmonics of Gaia baseline
        const_per_a_inds_list, const_per_m_inds_list, fmt =\
                                hlp_plot.period_lines(m_star, a_lim, m_lim,
                                                      a_lim, m_lim,
                                                      grid_num, 5, how='gaia')
        num_lines = len(const_per_a_inds_list)
        for i in range(num_lines):
            plt.plot(const_per_a_inds_list[i], const_per_m_inds_list[i], fmt, alpha=0.5, linewidth=3)
    
    
    ## Set plot limits ##
    # value2index(value, index_space, value_space)
    # x_min, x_max = value2index([a_lim], (0, grid_num-1), a_lim)
    # ax.set_xlim(a_lim)
    # ax.set_ylim(m_lim)
    # import pdb; pdb.set_trace()
    
    # a_min, m_min = a_lim[0], m_lim[0]
    # a_max, m_max = a_lim[1], m_lim[1]
    
    fig.tight_layout()
    save_dir = os.path.join(outdir, 'results/{}/'.format(star_name)) # Each star gets its own folder
    # Try to make directory. If it exists, just continue. Parallel code was bugging out here, so exist_ok is great.
    os.makedirs(save_dir, exist_ok = True)
    fig.savefig(save_dir + star_name + '_2d.png', dpi=400)
    
    plt.close()
    
    # # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # # twosig_inds contains the indices where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    # bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)
    
    return


def plot_1d(star_name, post_tot, a_lim, m_lim, 
            which=['cdf', 'pdf'], outdir=''):
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
        which (list of str): 'cdf' to plot/save CDFs, and 
                             'pdf to plot/save PDFs
        outdir (str): Path to save plot
    
    Returns:
        None (but plots and saves 1D posteriors)
        
    """
    bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)                
    
    title_size = 30
    label_size = 25
    tick_num = 6
    tick_size = 25
    
    sma_1d = post_tot.sum(axis=0)
    mass_1d = post_tot.sum(axis=1)
    
    grid_num = np.shape(post_tot)[0]
    
    for plot_type in which:
        fig, ax = plt.subplots(1,2, figsize=(12,8))
        
        if plot_type == 'cdf':
            plot_dist_a = np.cumsum(sma_1d)
            plot_dist_m = np.cumsum(mass_1d)
            
            max_ylim_a = 1
            max_ylim_m = 1
            
            ax[0].tick_params(axis='y', which='both', left=True, right=True, labelleft=True)
            ax[1].tick_params(axis='y', which='both', left=True, right=True, labelleft=True)
            
            save_name_suffix = '_cdf_1d.png'
            
        elif plot_type == 'pdf':
            plot_dist_a = sma_1d
            plot_dist_m = mass_1d
            
            max_ylim_a = np.max(sma_1d)
            max_ylim_m = np.max(mass_1d)
            
            ax[0].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
            ax[1].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
            
            save_name_suffix = '_pdf_1d.png'

        ax[0].set_ylim(0, max_ylim_a*1.05)
        ax[0].plot(range(grid_num + 1), np.insert(plot_dist_a, 0, 0))
        ax[0].set_title(f"Semi-major axis {plot_type.upper()}", size=title_size)
        ax[0].set_xlabel("Companion semi-major axis (AU)", size=label_size)

        ax[0].vlines(twosig_inds[0][0], 0, max_ylim_a, colors='r', linestyles='dashed')
        ax[0].vlines(twosig_inds[0][1], 0, max_ylim_a, colors='r', linestyles='dashed')
        ax[0].tick_params(axis='x', which='both', top=False, bottom=True)
        # ax[0].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)

        ax[1].set_ylim(0, max_ylim_m*1.05)
        ax[1].plot(range(grid_num + 1), np.insert(plot_dist_m, 0, 0))
        ax[1].set_title('Mass {}'.format(plot_type.upper()), size=title_size)
        ax[1].set_xlabel(r"Companion mass ($\mathrm{M_{Jup}}$)", size=label_size)

        ax[1].vlines(twosig_inds[1][0], 0, max_ylim_m, colors='r', linestyles='dashed')
        ax[1].vlines(twosig_inds[1][1], 0, max_ylim_m, colors='r', linestyles='dashed')
        ax[1].tick_params(axis='x', which='both', top=False, bottom=True)
        #ax[1].tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
        
        ###################################################
        ############ Tick positions and labels ################
        tick_num = 6
        min_exp = -4
        max_exp = 13
        n = max_exp-min_exp+1
        
        for i in range(2):
            lim = [a_lim, m_lim][i]
    
            ##
            tick_vals_long = np.logspace(min_exp, max_exp, n, base=10)

            # Exclude any positions that fall outside the grid domain
            tick_vals = tick_vals_long[(lim[0]<=tick_vals_long)\
                                          &(tick_vals_long<lim[1])]

            # Convert the labels to index positions. Note that the positions need not be integers, even though they correspond to "indices"
            tick_positions = hlp.value2index(tick_vals, (0, grid_num-1), lim)

            ax[i].set_xticks(tick_positions)
            ax[i].set_xticklabels([rf"$10^{{{int(np.log10(v))}}}$" for v in tick_vals], size=tick_size)
            
            ## Minor ticks ##
            minor_tick_list_of_arrs = [np.arange(2, 10) * 10.0**p for p in range(min_exp, max_exp)]
            minor_tick_vals_long = np.concatenate(minor_tick_list_of_arrs)
    
            # Exclude any positions that fall outside the grid domain
            minor_tick_vals = minor_tick_vals_long[(lim[0]<=minor_tick_vals_long)\
                                                      &(minor_tick_vals_long<lim[1])]
    
            minor_tick_positions = hlp.value2index(minor_tick_vals, (0, grid_num-1), lim)
    
            ax[i].set_xticks(minor_tick_positions, minor=True)
    
        save_dir = os.path.join(outdir, 'results/{}/'.format(star_name)) # Each star gets its own folder
        os.makedirs(save_dir, exist_ok = True)
    
        fig.tight_layout()
        fig.savefig(save_dir + star_name + save_name_suffix, dpi=400)
    
    return

