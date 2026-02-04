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

from astropy.time import Time


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
M_earth = 5.972167867791379e+27


def run(config_path, read_file_path=None, 
        plot=None, verbose=False):
    
    """
    Example API function to run orbit modeling and
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
        optional_params = ['num_points', 'min_a', 'max_a', 
                           'e_prior', 'a_m_prior', 'min_m', 'max_m',
                           'age_table', 'save', 'outdir']
        default_values = [int(1e6), 1, 64, 
                          'piecewise', 'cls', 1, 1e3, 4, ['proc'], '']

        num_points, min_a, max_a,\
        e_prior, a_m_prior, min_m, max_m,\
        age_table, save, outdir = driver.set_values(config_path, 
                                                    optional_params, 
                                                    default_values)
        
        num_points = int(num_points)
        ######################################
        ## Next load required params from config module
        cm = driver.load_module_from_file(config_path)
    
        star_name = cm.star_name
        m_star = cm.m_star
        d_star = cm.d_star
        grid_num = cm.grid_num
        
        grid_num = int(grid_num)
        ######################################
        
        ### General ###
        a_lim = (min_a, max_a)
        m_lim = (min_m, max_m)
    
        if verbose:
            print('Min sampling m is: ', min_m)
            print('Min sampling a is: ', min_a)
            
        ## Time the array calculations
        start_time = time.time()
        ##

        start_list_time = time.time()#######################################################
        a_list, m_list, per_list, e_list, i_list,\
        om_list, M_anom_0_list, a_inds, m_inds, log_a_m_prior = hlp.make_arrays(m_star, a_lim, m_lim,\
                                                                                grid_num, num_points,\
                                                                                e_prior, a_m_prior)
        
        end_list_time = time.time()#######################################################
                                                                 
        if verbose:
            print('made arrays')

        #######################################################################################
        ## RVs
        #######################################################################################
        start_rv_time = time.time()#######################################################
        if cm.run_rv:
            gammadot = cm.gammadot
            gammadot_err = cm.gammadot_err
            gammaddot = cm.gammaddot
            gammaddot_err = cm.gammaddot_err
            rv_epoch = cm.rv_epoch
            
            if (gammaddot is None) or (gammaddot_err is None):
                gammaddot, gammaddot_err = 0, 1e8
                
            # rv_list is a 1D array of log-likelihoods for the RV data
            # post_rv is the result of exponentiating rv_list and reshaping it into a 2D array
            rv_list = hlp_rv.rv_list(a_list, m_list, e_list, i_list, om_list, M_anom_0_list,
                                     per_list, cm.m_star, rv_epoch,
                                     gammadot, gammadot_err, gammaddot, gammaddot_err)

            post_rv = hlp.post_single(rv_list, log_a_m_prior, a_inds, m_inds, grid_num)
            
        else:
            rv_list = np.zeros(num_points) # 1D arrays contain log-likelihoods
            post_rv = np.ones((grid_num, grid_num)) # 2D arrays contain likelihoods (actually posteriors due to binning)
        end_rv_time = time.time()#######################################################
        #######################################################################################
        ## Astrometry
        #######################################################################################
        start_astro_time = time.time()#######################################################
        if cm.run_astro:
            delta_mu = cm.delta_mu
            delta_mu_err = cm.delta_mu_err
            hip_id = cm.hip_id
            gaia_id = cm.gaia_id
        
            # If delta_mu is not provided directly, use provided name
            if any([val is None for val in [delta_mu, delta_mu_err]]):
                delta_mu, delta_mu_err = hlp_astro.HGCA_retrieval(hip_id, gaia_id)
                
            # astro_list is a 1D array of log-likelihoods for the astrometry data
            # post_astro is the result of exponentiating astro_list and reshaping it into a 2D array
            astro_list = hlp_astro.astro_list(a_list, m_list, e_list, i_list, 
                                              om_list, M_anom_0_list, per_list,
                                              m_star, d_star, delta_mu, delta_mu_err)                     

            post_astro = np.array(hlp.post_single(astro_list, log_a_m_prior, a_inds, m_inds, grid_num))

        else:
            astro_list = np.zeros(num_points) # 1D arrays contain log-likelihoods
            post_astro = np.ones((grid_num, grid_num)) # 2D arrays contain likelihoods
        end_astro_time = time.time()#######################################################
        #######################################################################################
        ## Imaging
        #######################################################################################
        start_imag_time = time.time()#######################################################
        if cm.run_imag:
            imag_calc = driver.set_values(config_path, ['imag_calc'], ['exact'])
            vmag = cm.vmag
            imag_wavelength = cm.imag_wavelength
            contrast_str = cm.contrast_str
        
        
            if imag_calc == 'exact':
                imag_epoch = cm.imag_epoch
                
                # imag_list is a 1D array of log-likelihoods for the imaging data
                # post_imag is the result of exponentiating imag_list and reshaping it into a 2D array
                imag_list = hlp_imag.imag_list(a_list, m_list, e_list, i_list, om_list, 
                                               M_anom_0_list, per_list, m_star, 
                                               d_star, vmag, imag_wavelength, age_table,
                                               imag_epoch, contrast_str)
                                               
                #post_imag= hlp.post_single(imag_list, log_a_m_prior, a_inds, m_inds, grid_num)
                ## New change: save and plot the *approx* imag array even when calculating imag exactly (but still calculate post_tot exactly)
                post_imag = hlp_imag.imag_array(d_star, vmag, imag_wavelength, age_table,
                                                contrast_str, a_lim, m_lim, grid_num)
    
            elif imag_calc == 'approx':
                imag_list = np.zeros(num_points) # Dummy list to pass to tot_list() function
                post_imag = hlp_imag.imag_array(d_star, vmag, imag_wavelength, age_table,
                                                contrast_str, a_lim, m_lim, grid_num)
            
            else:
                raise Exception("api_run.run: 'imag_calc' must be either 'exact' or 'approx'")
    
        else:
            imag_list = np.zeros(num_points) # 1D arrays contain log-likelihoods
            post_imag = np.ones((grid_num, grid_num)) # 2D arrays contain likelihoods
            vmag=None
            imag_wavelength=None
            contrast_str=None
            imag_calc=None
        end_imag_time = time.time()#######################################################
        ##################################################################################
        ## Total
        ###################################################################################
        start_tot_time = time.time()#######################################################
        
        # tot_list is a 1D array of log-likelihoods for all data. No priors incorporated yet.
        tot_list = np.array(hlp.tot_list(rv_list, astro_list, imag_list, num_points))
        
        # If imag=='exact', then post_tot is the result of exponentiating imag_list and reshaping it into a 2D array. Otherwise there is a special treatment for the imag_calc=="approx" case.
        if cm.run_imag and imag_calc=='exact':
            post_tot = hlp.post_single(tot_list, log_a_m_prior, a_inds, m_inds, grid_num)
        
        else:
            # Need to multiply by post_imag in approx imaging case
            post_tot = hlp.post_tot_approx_imag(tot_list, post_imag, log_a_m_prior, a_inds, m_inds, grid_num)
            
        end_tot_time = time.time()#######################################################
        #######################################################################################
        #######################################################################################

        ##
        end_time = time.time()
        ##
        tot = end_time-start_time
        
        
        if verbose:
            
            print('{:.0e} points ran for {} in {:.2f} seconds.'.format(num_points, star_name, end_time-start_time))
            print('')
            print('Fractions of total time:')
            print('Arrays: {:.2f}'.format((end_list_time-start_list_time)/tot))
            print('RVs: {:.2f}'.format((end_rv_time-start_rv_time)/tot))
            print('Astro: {:.2f}'.format((end_astro_time-start_astro_time)/tot))
            print('Imag: {:.2f}'.format((end_imag_time-start_imag_time)/tot))
            print('Creating post_tot: {:.2f}'.format((end_tot_time-start_tot_time)/tot))
        

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
                        tot_list, rv_list, astro_list, imag_data,
                        vmag, imag_wavelength, contrast_str, age_table,
                        log_a_m_prior, a_list, m_list,
                        a_lim, m_lim, imag_calc=imag_calc, outdir=outdir, 
                        verbose=False)

    # If read_file_path is NOT None, load in existing data:
    else:
        # Load grid_num from config file
        cm = driver.load_module_from_file(config_path)
        
        star_name, m_star, d_star,\
        run_rv, run_astro, run_imag,\
        post_tot, post_rv, post_astro, post_imag,\
        a_lim, m_lim = ls.load(read_file_path, grid_num=cm.grid_num, verbose=verbose)
        
        
    if plot==True:
        
        # Check if scatter_plot and outdir are provided in config. Otherwise set to defaults
        scatter_plot, outdir = driver.set_values(config_path, ['scatter_plot', 'outdir'], [None, ''])
        
        plotter.joint_plot(star_name, m_star, d_star,
                           run_rv, run_astro, run_imag,
                           post_tot, post_rv, post_astro, post_imag, 
                           a_lim, m_lim,
                           scatter_plot=scatter_plot, 
                           period_lines=False, outdir='', verbose=verbose)
        plotter.plot_1d(star_name, post_tot, a_lim, m_lim, outdir='')
    
        
    # Print out the 2-sigma boundaries (bounds) for the joint posterior
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained within the interval twosig_inds[i]
    if verbose:
        # bounds is the final answer: [range of 2σ a, range of 2σ m].
        # twosig_inds contains the indices corresponding to bounds. That is, where the CDF reaches the upper and lower values associated with the 95% confidence interval.
        bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)
        print("a_lim = ", bounds[0])
        print("m_lim = ", bounds[1])
    return


if __name__ == "__main__":

    config_path = 'test_config_files/test1.py'
    read_file_path = None#'results/test1/test1_raw.h5'


    plot=True
    verbose=True
    first_time = time.time()
    run(config_path, read_file_path,
        plot=plot, verbose=verbose)
    second_time = time.time()
    print("TIME {:.2f}".format(second_time-first_time))   
    
    
    
    
    
    
    
    
    
    
    

