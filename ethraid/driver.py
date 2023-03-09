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
        contrast_curve (dataframe or dict, 
                        columns of 'ang_sep' (arcseconds) and 'delta_mag' (mag)): 
                        Ordered pairs of angular separation and Δmag.
    
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
    scatter_plot = args.scatter_plot
    num_points = args.num_points
    grid_num = args.grid_num
    save = args.save
    plot = args.plot
    read_file_path = args.read_file_path
    outdir = args.outdir
    

    # If no data to read in, calculate new arrays
    if read_file_path is None:
        
        min_per = rv_baseline*0.7
        min_m = 1

        m_star_Ms = m_star * M_jup/M_sun
        # Finally, the minimum semi-major axis is the one where period is smallest and companion mass is smallest too. If companion mass were larger at the same period, the companion would have to be farther away. Same for larger period at fixed mass.
        min_a = rv.utils.semi_major_axis(min_per, (m_star_Ms + min_m*(M_jup/M_sun)))
        

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
            print('No astrometry data provided. Bounds will be based on RVs only.')
    

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

        if save==True:
            no_astro = True if (delta_mu is None or delta_mu_err is None) else False
                
            ls.save(star_name, m_star, d_star, rv_list, astro_list, no_astro, post_imag, a_list, m_list,
                    a_lim, m_lim, outdir=outdir)
    
    # Otherwise, load in existing data:
    else:
        star_name, m_star, d_star,\
        post_tot, post_rv, post_astro, post_imag,\
        a_lim, m_lim = ls.load(read_file_path, grid_num)
        
        
    if plot==True:
        plotter.joint_plot(star_name, m_star, d_star, vmag, 
                           post_tot, post_rv, post_astro, post_imag, 
                           grid_num, a_lim, m_lim, scatter_plot=scatter_plot, 
                           period_lines = False, outdir=outdir)
    
    return
    

if __name__ == "__main__":
    
    run(*sp.params_191939_old, num_points=1e6, grid_num=100, plot=True, read_file_path=None, outdir='')
    #'results/post_arrays/191939_old_1e8.h5')
    #'results/post_arrays/12572.h5'
    # run(*sp.params_synth, num_points=1e6, grid_num=100, save=False, plot=True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

