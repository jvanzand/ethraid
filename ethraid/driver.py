import astropy.constants as c
import numpy as np
import time
import radvel as rv

import os
import sys

from ethraid import plotter
from ethraid import system_params as sp
from ethraid import load_save as ls
from ethraid import helper_functions_imaging as hlp_imag

#########################
from ethraid.compiled import helper_functions_general as hlp
from ethraid.compiled import helper_functions_rv as hlp_rv
from ethraid.compiled import helper_functions_astro as hlp_astro
#########################


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
M_earth = 5.972167867791379e+27

def run(args):
        
    """
    Primary function to run code.
        
    Arguments:
        args: Command line arguments, including:
    
        star_name (str): Name of host star
        m_star (Jupiter masses): Mass of host star
        d_star (AU): Distance to host star
    
        gammadot (m/s/day): RV trend term
        gammadot_err (m/s/day): Error on gammadot
        gammaddot (m/s/day/day): RV curvature term
        gammaddot_err (m/s/day/day): Error on gammaddot
        rv_baseline (days): Time star has been observed
        rv_epoch (BJD): Date at which model gammadot and 
                        gammaddot will be evaluated. 
                        Typically ~halfway through the baseline
    
        delta_mu (milli-arcseconds/year): Magnitude of the change
                        in astrometric proper motion between 
                        Hipparcos and Gaia epochs.
        delta_mu_err: Error on delta_mu
    
        vmag (mag): Apparent V-band magnitude of host star
        imag_wavelength (μm): Wavelength of imaging observations
        contrast_str (dataframe or dict, 
                      columns of 'ang_sep' (arcseconds) 
                      and 'delta_mag' (mag)
                      ): Ordered pairs of angular separation and Δmag.
        
        num_points (int): Number of orbits to simulate (usually 1e6 - 1e8)
        grid_num (int): Dimensions of 2D posterior. Determines "resolution"
        
        save (bool): Whether to save raw probability arrays.
                     Processed arrays will be saved by default.
        outdir (str): Path to save outputs to
        verbose (bool): Optionally print out extra information
    
    Returns:
        None
    """

    star_name = args.star_name
    m_star = args.m_star
    d_star = args.d_star
    gammadot = args.gammadot
    gammadot_err = args.gammadot_err
    gammaddot = args.gammaddot
    gammaddot_err = args.gammaddot_err
    rv_baseline = args.rv_baseline
    rv_epoch = args.rv_epoch
    delta_mu = args.delta_mu
    delta_mu_err = args.delta_mu_err
    vmag = args.vmag
    imag_wavelength = args.imag_wavelength
    contrast_str = args.contrast_str
    num_points = args.num_points
    grid_num = args.grid_num
    save = args.save
    outdir = args.outdir
    verbose = args.verbose

        
    min_per = rv_baseline*2
    min_m = hlp.min_mass(gammadot, gammaddot, rv_baseline, min_per, m_star)
    # Finally, the minimum semi-major axis is the one where period is smallest and companion mass is smallest too. If companion mass were larger at the same period, the companion would have to be farther away. Same for larger period at fixed mass.
    # Radvel requires solar masses.
    min_a = rv.utils.semi_major_axis(min_per, ((m_star + min_m)*(M_jup/M_sun)))
    

    print('Min sampling m is: ', min_m)
    print('Min sampling a is: ', min_a)
    
    
    ### General ###
    # Arbitrary upper limits
    max_a = 1e2
    max_m = 5e2
    
    a_lim = (min_a, max_a)
    m_lim = (min_m, max_m)

    num_points = int(num_points)
    grid_num = int(grid_num)

    a_list, m_list, per_list, e_list, i_list,\
    om_list, M_anom_0_list, a_inds, m_inds = hlp.make_arrays(m_star, a_lim, m_lim,\
                                                            grid_num, num_points)

    print('made arrays')
    ##
    start_time = time.time()
    ##
    
    ## Start with the imaging posterior. This rules out any companions massive enough to be visible in imaging data.
    post_imag = hlp_imag.imag_array(d_star, vmag, imag_wavelength, 
                                    contrast_str, a_lim, m_lim, grid_num)
    
    ## Now the astrometry posterior.
    # Some targets aren't in the Hip/Gaia catalog, so we can't make the astrometry posterior for them.
    try:
        astro_list = hlp_astro.astro_list(a_list, m_list, e_list, i_list, 
                                          om_list, M_anom_0_list, per_list,
                                          m_star, d_star, delta_mu, delta_mu_err)                     
                             
        post_astro = np.array(hlp.post_single(astro_list, a_inds, m_inds, grid_num))


    except Exception as err:
        astro_list = np.ones(num_points)
        post_astro = np.zeros((grid_num, grid_num))
        
        if verbose:
            print('driver.run: No astrometry data provided. Bounds will be based on RVs only.')


    ## Last we calculate the RV posterior
    rv_list = hlp_rv.rv_list(a_list, m_list, e_list, i_list, om_list, 
                             M_anom_0_list, per_list, m_star, rv_epoch,
                             gammadot, gammadot_err, gammaddot, gammaddot_err)
                            
    post_rv = hlp.post_single(rv_list, a_inds, m_inds, grid_num)
    
    post_tot = hlp.post_tot(rv_list, astro_list, post_imag, grid_num, a_inds, m_inds)
    
    ##
    end_time = time.time()
    ##
    print('{:.0e} points ran for {} in {:.2f} seconds.'.format(num_points, star_name, end_time-start_time))

    no_astro = True if None in [delta_mu, delta_mu_err] else False
    no_imag = True if None in [vmag, imag_wavelength, contrast_str] else False
    
    if 'proc' in save:
        ls.save_processed(star_name, m_star, d_star, post_tot,
                          post_rv, post_astro, post_imag, 
                          no_astro, no_imag,
                          a_lim, m_lim, outdir=outdir)
    if 'raw' in save:
        ls.save_raw(star_name, m_star, d_star, 
                    rv_list, astro_list, post_imag, 
                    no_astro, no_imag,
                    vmag, imag_wavelength, contrast_str,
                    a_list, m_list, a_lim, m_lim, 
                    outdir=outdir, verbose=verbose)
        
    return
    
    
def plot(args):
    """
    Plot the content of loaded or calculated arrays
    
    Arguments:
        args: Command line arguments, including:
        
        read_file_path (str): Path to file containing data to plot
        type (str): Either '2d' or '1d'
        grid_num (int): Dimensions of 2D probability array
                        This argument need not be provided IF
                        read_file_path points to processed arrays
                        rather than raw arrays.
        outdir (str): Path to which plots should be saved
    
    Returns:
        None
    """
    
    star_name, m_star, d_star,\
    post_tot, post_rv, post_astro, post_imag,\
    grid_num, a_lim, m_lim = ls.load(args.read_file_path, args.grid_num, args.verbose)
    
    if "2d" in args.type:
        plotter.joint_plot(star_name, m_star, d_star, 
                           post_tot, post_rv, post_astro, post_imag, 
                           grid_num, a_lim, m_lim, 
                           scatter_plot=args.scatter_plot, period_lines=False,
                           outdir=args.outdir, verbose=args.verbose)
    
    if "1d" in args.type:
        plotter.plot_1d(star_name, post_tot, a_lim, m_lim, outdir=args.outdir)
                       
    return
    
    
def less(args):
    """
    Load in saved arrays and calculate 1D a and m bounds.
    Print out and return these bounds.
    
    Arguments:
        args: Command line arguments, including:

        read_file_path (str): Path to saved outputs.

    Returns:
        bounds (list of tuples): [(a1, a2),(m1, m2)] giving 95%
                                 confidence intervals
    """ 
    

    star_name, m_star, d_star,\
    post_tot, post_rv, post_astro, post_imag,\
    grid_num, a_lim, m_lim = ls.load(args.read_file_path, args.grid_num, args.verbose)
    
    # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # twosig_inds contains the indices corresponding to bounds. That is, where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)
        
    # Print out the 2-sigma boundaries (bounds) for the joint posterior
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained within the interval twosig_inds[i]
    print('a_lim = ', bounds[0], ' AU')
    print('m_lim = ', bounds[1], ' M_J')
    
    return bounds
    
    
    
    
    
    
    
    
    

