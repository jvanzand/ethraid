## This is a version of driver.py that can be run directly as a module.

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

def run(star_name=None, m_star=None, d_star=None, 
        gammadot=None, gammadot_err=None, 
        gammaddot=0, gammaddot_err=1e8, 
        rv_baseline=None, rv_epoch=None, 
        delta_mu=None, delta_mu_err=None,
        vmag=None, imag_wavelength=None, contrast_str=None, 
        scatter_plot=None, num_points=1e6, grid_num=100, 
        plot=True, read_file_path=None, save=False, 
        outdir='', verbose=False):
        
    """
    Example API function to run orbit modelling and
    save and plot results.
        
    Arguments:
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
        contrast_curve (dataframe or dict, 
                        columns of 'ang_sep' (arcseconds) and 'delta_mag' (mag)): 
                        Ordered pairs of angular separation and Δmag.
        
        scatter_plot (2-tuple of floats): (a, m) of known companion to plot on 2D
                                          posterior
        num_points (int): Number of orbits to simulate (usually 1e6 - 1e8)
        grid_num (int): Dimensions of 2D posterior. Determines "resolution"
        
        read_file_path (str): Path to saved outputs. Providing this 
                              will circumvent core probability 
                              calculations.
        save (bool): Whether to save raw probability arrays.
                     Processed arrays will be saved by default.
        outdir (str): Path to save outputs to
        verbose (bool): Optionally print out extra information
    
    Returns:
        None
    """
    if read_file_path is None and None in [star_name, m_star, d_star, gammadot, gammadot_err]:
        raise Exception("Either read_file_path or all of \n"
                        "          {star_name, m_star, d_star, gammadot, gammadot_err, \n"
                        "          rv_baseline, rv_epoch} must be provided.")

    # If no data to read in, calculate new arrays
    if read_file_path is None:
        
        min_per = rv_baseline*2
        min_m = hlp.min_mass(gammadot, gammaddot, rv_baseline, min_per, m_star)
        # Finally, the minimum semi-major axis is the one where period is smallest and companion mass is smallest too. If companion mass were larger at the same period, the companion would have to be farther away. Same for larger period at fixed mass.
        min_a = rv.utils.semi_major_axis(min_per, ((m_star + min_m)*(M_jup/M_sun)))

        print('Min sampling m is: ', min_m)
        print('Min sampling a is: ', min_a)

        max_a = 1e2
        max_m = 1e3
        
        ### General ###
        # Arbitrary upper limits
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
        post_imag = hlp_imag.imag_array(d_star, vmag, imag_wavelength, contrast_str, a_lim, m_lim, grid_num)
        
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
                print('api_run.run: No astrometry data provided. Bounds will be based on RVs only.')
    

        ## Last we calculate the RV posterior
        rv_list = hlp_rv.rv_list(a_list, m_list, e_list, i_list, om_list, M_anom_0_list,
                                per_list, m_star, rv_epoch,
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
    
    # Otherwise, load in existing data:
    else:
        star_name, m_star, d_star,\
        post_tot, post_rv, post_astro, post_imag,\
        grid_num, a_lim, m_lim = ls.load(read_file_path, grid_num, verbose)
        
        
    if plot==True:
        plotter.joint_plot(star_name, m_star, d_star, 
                           post_tot, post_rv, post_astro, post_imag, 
                           grid_num, a_lim, m_lim, scatter_plot=scatter_plot, 
                           period_lines = False, outdir=outdir, verbose=verbose)
        plotter.plot_1d(star_name, post_tot, a_lim, m_lim, outdir=outdir)
    
    # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # twosig_inds contains the indices corresponding to bounds. That is, where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)
        
    # Print out the 2-sigma boundaries (bounds) for the joint posterior
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained within the interval twosig_inds[i]
    print('a_lim = ', bounds[0], ' AU')
    print('m_lim = ', bounds[1], ' M_J')
    
    return
    

if __name__ == "__main__":
    
    rfp = 'results/191939/191939_raw.h5'
    run(*sp.params_191939, num_points=1e6, grid_num=100, plot=True, read_file_path=None, 
        save=['proc', 'raw'], outdir='', verbose=True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

