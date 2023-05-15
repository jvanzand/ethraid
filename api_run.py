## This is a version of driver.py that can be run directly as a module.
import astropy.constants as c
import numpy as np
import time

import os
import sys

from ethraid import plotter
from ethraid import load_save as ls
from ethraid import driver

#########################
from ethraid.compiled import helper_functions_general as hlp
from ethraid.compiled import helper_functions_rv as hlp_rv
from ethraid.compiled import helper_functions_astro as hlp_astro
from ethraid.compiled import helper_functions_imaging as hlp_imag
#########################


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
M_earth = 5.972167867791379e+27


def run(config_path=None, read_file_path=None, 
        plot=None, verbose=False):
    
    """
    Example API function to run orbit modelling and
    save and plot results.
        
    Arguments:
        config_path (str): Path to configuration file.
                           Will only be read if read_file_path
                           is None.
        read_file_path (str): Path to saved results (.h5 file)
                              If provided, will override params
                              given in config_path.
        
        The arguments below are only needed if read_file_path is not None.
        If you are instead using a config file, define them there (those
        provided in the config file will supersede parameters provided 
        directly).
        
        plot (bool): Plot results?
        verbose (bool): Verbose output?
    
    Returns:
        None
    """

    # If no data to read in, calculate new arrays
    if read_file_path is None:

        ## First try loading optional params. If they can't be loaded, use defaults
        optional_params = ['num_points', 'grid_num', 'save', 'min_a', 'max_a', 'min_m', 'max_m']
        default_values = [int(1e6), int(1e2), ['proc'], 1, 1e2, 1, 1e3]

        num_points, grid_num, save,\
        min_a, max_a, min_m, max_m = driver.set_values(config_path, 
                                                       optional_params, 
                                                       default_values)
        
        num_points = int(num_points)
        grid_num = int(grid_num)
        ######################################
        ## Next load in required params
        cm = driver.load_module_from_file(config_path)
    
        star_name = cm.star_name
        m_star = cm.m_star
        d_star = cm.d_star
        ######################################
        
        ### General ###
        a_lim = (min_a, max_a)
        m_lim = (min_m, max_m)
    
        if verbose:
            print('Min sampling m is: ', min_m)
            print('Min sampling a is: ', min_a)
            
        ## Time array calculations
        start_time = time.time()
        ##
        
        start_time_arrays = time.time()
        a_list, m_list, per_list, e_list, i_list,\
        om_list, M_anom_0_list, a_inds, m_inds = hlp.make_arrays(cm.m_star, a_lim, m_lim,\
                                                                 grid_num, num_points)
        end_time_arrays = time.time()
                                                                 
        if verbose:
            print('made arrays')
            
        start_time_rv = time.time()
        #######################################################################################
        ## RVs
        #######################################################################################
        if cm.run_rv:
            gammadot = cm.gammadot
            gammadot_err = cm.gammadot_err
            gammaddot = cm.gammaddot
            gammaddot_err = cm.gammaddot_err
            rv_epoch = cm.rv_epoch
        
            rv_list = hlp_rv.rv_list(a_list, m_list, e_list, i_list, om_list, M_anom_0_list,
                                    per_list, cm.m_star, rv_epoch,
                                    gammadot, gammadot_err, gammaddot, gammaddot_err)
            post_rv = hlp.post_single(rv_list, a_inds, m_inds, grid_num)
    
        else:
            rv_list = np.ones(num_points)
            post_rv = np.ones((grid_num, grid_num))
        end_time_rv = time.time()
        #######################################################################################
        ## Astrometry
        #######################################################################################
        start_time_astro = time.time()
        if cm.run_astro:
            delta_mu = cm.delta_mu
            delta_mu_err = cm.delta_mu_err
            hip_id = cm.hip_id
            gaia_id = cm.gaia_id
        
            # If delta_mu is not provided directly, use provided name
            if any([val is None for val in [delta_mu, delta_mu_err]]):
                delta_mu, delta_mu_err = hlp_astro.HGCA_retrieval(hip_id, gaia_id)
    
            astro_list = hlp_astro.astro_list(a_list, m_list, e_list, i_list, 
                                              om_list, M_anom_0_list, per_list,
                                              m_star, d_star, delta_mu, delta_mu_err)                     

            post_astro = np.array(hlp.post_single(astro_list, a_inds, m_inds, grid_num))

        else:
            astro_list = np.ones(num_points)
            post_astro = np.ones((grid_num, grid_num))
        end_time_astro = time.time()
        #######################################################################################
        ## Imaging
        #######################################################################################
        start_time_imag = time.time()
        if cm.run_imag:
            imag_calc = driver.set_values(config_path, ['imag_calc'], ['exact'])
            vmag = cm.vmag
            imag_wavelength = cm.imag_wavelength
            contrast_str = cm.contrast_str
        
        
            if imag_calc == 'exact':
                imag_epoch = cm.imag_epoch
                
                imag_list = hlp_imag.imag_list(a_list, m_list, e_list, i_list, om_list, 
                                               M_anom_0_list, per_list, m_star, 
                                               d_star, vmag, imag_wavelength, 
                                               imag_epoch, contrast_str)
                post_imag= hlp.post_single(imag_list, a_inds, m_inds, grid_num)
    
            elif imag_calc == 'approx':
                post_imag = hlp_imag.imag_array(d_star, vmag, imag_wavelength, 
                                                contrast_str, a_lim, m_lim, grid_num)
    
        else:
            imag_list = np.ones(num_points)
            post_imag = np.ones((grid_num, grid_num))
            vmag=None
            imag_wavelength=None
            contrast_str=None
            imag_calc=None
        end_time_imag = time.time()   
        #######################################################################################
        ## Total
        #######################################################################################
        start_time_post_tot = time.time()
        if cm.run_imag and imag_calc=='exact':
            post_tot = hlp.post_tot(rv_list, astro_list, imag_list, grid_num, a_inds, m_inds)
        
        else:
            post_tot = hlp.post_tot_simplified(rv_list, astro_list, post_imag, grid_num, a_inds, m_inds)
        end_time_post_tot = time.time()
        #######################################################################################
        #######################################################################################
        
        ##
        end_time = time.time()
        ##
        
        if verbose:
            # print('Arrays {:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time_arrays-start_time_arrays))
            #
            # print('RV {:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time_rv-start_time_rv))
            # print('Astro {:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time_astro-start_time_astro))
            # print('Imag {:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time_imag-start_time_imag))
            
            
            
            print('{:.0e} points ran for {} in {:.2f} seconds.'.format(num_points, star_name, end_time-start_time))
        
        #################################################
        time_arrays = end_time_arrays-start_time_arrays
        time_rv = end_time_rv-start_time_rv
        time_astro = end_time_astro-start_time_astro
        time_imag = end_time_imag-start_time_imag
        time_post_tot = end_time_post_tot-start_time_post_tot
        
        print('IMAGING TOOK {} SECONDS'.format(time_imag))
        
        time_tot = time_arrays + time_rv + time_astro + time_imag + time_post_tot
        t_types = ['Array', 'RV', 'Astro', 'Imag', 'Post_tot']
        print('')
        for i, t in enumerate([time_arrays, time_rv, time_astro, time_imag, time_post_tot]):
            print("{}'s fraction of total time: {:.2f}".format(t_types[i], t/time_tot) )
        ##################################################
        run_rv = cm.run_rv
        run_astro = cm.run_astro
        run_imag = cm.run_imag
        
        
        # Check if scatter_plot and outdir are provided in config. Otherwise set to defaults
        scatter_plot, outdir = driver.set_values(config_path, ['scatter_plot', 'outdir'], [None, ''])
        
        if 'proc' in save:
            ls.save_processed(star_name, m_star, d_star,
                              run_rv, run_astro, run_imag, 
                              post_tot, post_rv, post_astro, post_imag,
                              a_lim, m_lim, outdir=outdir)

        if 'raw' in save:
            if imag_calc=='approx':
                imag_data = post_imag
            else: # If cm.imag_calc='exact' or anything else
                imag_data = imag_list
            
            ls.save_raw(star_name, m_star, d_star, 
                        run_rv, run_astro, run_imag,
                        rv_list, astro_list, imag_data,
                        vmag, imag_wavelength, contrast_str,
                        a_list, m_list, a_lim, m_lim, 
                        imag_calc=imag_calc, outdir=outdir, 
                        verbose=False)

    # If read_file_path is NOT None, load in existing data:
    else:
        star_name, m_star, d_star,\
        run_rv, run_astro, run_imag,\
        post_tot, post_rv, post_astro, post_imag,\
        grid_num, a_lim, m_lim = ls.load(read_file_path, grid_num, verbose)
        
        
    if plot==True:
        plotter.joint_plot(star_name, m_star, d_star,
                           run_rv, run_astro, run_imag,
                           post_tot, post_rv, post_astro, post_imag, 
                           grid_num, a_lim, m_lim,
                           scatter_plot=scatter_plot, 
                           period_lines = False, outdir=outdir, verbose=verbose)
        plotter.plot_1d(star_name, post_tot, a_lim, m_lim, outdir=outdir)
    
    # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # twosig_inds contains the indices corresponding to bounds. That is, where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)
        
    # Print out the 2-sigma boundaries (bounds) for the joint posterior
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained within the interval twosig_inds[i]
    if verbose:
        print('a_lim = ', bounds[0], ' AU')
        print('m_lim = ', bounds[1], ' M_J')
    
    return time_arrays, time_rv, time_astro, time_imag, time_post_tot, time_tot


if __name__ == "__main__":
    
    config_path = 'ethraid/example_config_files/config_191939.py'
    read_file_path = None
    
    plot=True
    verbose = True
    
    run(config_path, read_file_path,
        plot=plot, verbose=verbose)
    
    
    
    
    
    
    
    
    
    
    
    

