import astropy.constants as c
import numpy as np
import time

import os
import sys

from ethraid import plotter
from ethraid import load_save as ls

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


def run(args):
        
    """
    Primary function to run orbit fitting code.
        
    Arguments:
        args: Command line arguments, especially config,
              a config file which contains the following:
    
        star_name (str): Name of host star
        m_star (Jupiter masses): Mass of host star
        d_star (AU): Distance to host star
    
        run_rv (bool): Use RV data in calculation?
        run_astro (bool): Use astrometry data in calculation?
        run_imag (bool): Use imaging data in calculation?
    
        gammadot (m/s/day): RV trend term
        gammadot_err (m/s/day): Error on gammadot
        gammaddot (m/s/day/day): RV curvature term
        gammaddot_err (m/s/day/day): Error on gammaddot
        rv_epoch (BJD): Date at which model gammadot and 
                        gammaddot will be evaluated. 
                        Typically ~halfway through the baseline
    
        delta_mu (milli-arcseconds/year): Magnitude of the change
                        in astrometric proper motion between 
                        Hipparcos and Gaia epochs.
        delta_mu_err: Error on delta_mu
        hip_id (string): Hipparcos identifier to retrieve astrometry data;
                         alternative to providing delta_mu/delta_mu_err.
                         If hip_id and gaia_id are both provided, they
                         must correspond to the same target.
        gaia_id (string): Gaia identifier to retrieve astrometry data;
                          alternative to providing delta_mu/delta_mu_err.
                          If hip_id and gaia_id are both provided, they
                          must correspond to the same target.
    
        vmag (mag): Apparent V-band magnitude of host star
        imag_wavelength (μm): Wavelength of imaging observations
        age_table (int): Integer 1-5, indicating which BD cooling model to use
                         based on age of system.
                         1-->0.1 Gyr, 2-->0.5 Gyr, 3-->1 Gyr, 4-->5 Gyr, 5-->10 Gyr
        contrast_str (dataframe or dict, 
                      columns of 'ang_sep' (arcseconds) 
                      and 'delta_mag' (mag)
                      ): Ordered pairs of angular separation and Δmag.
        
        num_points (int): Number of orbits to simulate (usually 1e6 - 1e8)
        grid_num (int): Dimensions of 2D posterior. Determines "resolution"
        
        save (list): List containing 'raw' and/or 'proc'
                     to indicate whether raw or processed arrays
                     should be saved. With the CLI, only processed 
                     arrays will be saved by default.
        outdir (str): Path to save outputs to
        scatter_plot (list): Optional list of (sep, mass) tuples to scatter plot 
                             the parameters of 1 or more companions. Sma in AU,
                             mass in M_jup.
    
    Returns:
        None
    """
    
    config_path = args.config
    
    ## First try loading optional params. If they can't be loaded, use defaults
    optional_params = ['num_points', 'min_a', 'max_a', 
                       'e_dist', 'min_m', 'max_m', 'save', 'outdir']
    default_values = [1e6, 1, 1e2, 
                      'piecewise', 1, 1e3, ['proc'], '']

    num_points, min_a, max_a,\
    e_dist, min_m, max_m, save, outdir = set_values(config_path, 
                                            optional_params, 
                                            default_values)
    
    num_points = int(num_points)
    ######################################
    ## Next load required params from config module
    cm = load_module_from_file(config_path)

    star_name = cm.star_name
    m_star = cm.m_star
    d_star = cm.d_star
    grid_num = cm.grid_num
    verbose = args.verbose
    
    
    run_rv = cm.run_rv
    run_astro = cm.run_astro
    run_imag = cm.run_imag
    
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

    a_list, m_list, per_list, e_list, i_list,\
    om_list, M_anom_0_list, a_inds, m_inds, prior = hlp.make_arrays(cm.m_star, a_lim, m_lim,\
                                                                    grid_num, num_points, e_dist)

    if verbose:
        print('made arrays')
    
    #######################################################################################
    ## RVs
    #######################################################################################
    if cm.run_rv:
        gammadot = cm.gammadot
        gammadot_err = cm.gammadot_err
        gammaddot = cm.gammaddot
        gammaddot_err = cm.gammaddot_err
        rv_epoch = cm.rv_epoch
        
        if (gammaddot is None) or (gammaddot_err is None):
            gammaddot, gammaddot_err = 0, 1e8
        
        rv_list = hlp_rv.rv_list(a_list, m_list, e_list, i_list, om_list, M_anom_0_list,
                                per_list, cm.m_star, rv_epoch,
                                gammadot, gammadot_err, gammaddot, gammaddot_err)
        post_rv = hlp.post_single(rv_list, a_inds, m_inds, grid_num)
    
    # If run_rv is False, populate arrays with 1s
    else:
        rv_list = np.zeros(num_points)
        post_rv = np.ones((grid_num, grid_num))
    
    #######################################################################################
    ## Astrometry
    #######################################################################################
    if run_astro:
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
        astro_list = np.zeros(num_points)
        post_astro = np.ones((grid_num, grid_num))
    
    #######################################################################################
    ## Imaging
    #######################################################################################
    if run_imag:
        imag_calc = set_values(config_path, ['imag_calc'], ['exact'])
        vmag = cm.vmag
        imag_wavelength = cm.imag_wavelength
        age_table = cm.age_table
        contrast_str = cm.contrast_str
        
        
        if imag_calc == 'exact':
            imag_epoch = cm.imag_epoch
            
            imag_list = hlp_imag.imag_list(a_list, m_list, e_list, i_list, om_list, 
                                           M_anom_0_list, per_list, m_star, 
                                           d_star, vmag, imag_wavelength, age_table, 
                                           imag_epoch, contrast_str)
            post_imag= hlp.post_single(imag_list, a_inds, m_inds, grid_num)
    
        elif imag_calc == 'approx':
            imag_list = np.zeros(num_points) # Dummy list to pass to tot_list() function
            post_imag = hlp_imag.imag_array(d_star, vmag, imag_wavelength, age_table, 
                                            contrast_str, a_lim, m_lim, grid_num)
        
        else:
            raise Exception("driver.run: 'imag_calc' must be either 'exact' or 'approx'")
    
    else:
        imag_list = np.zeros(num_points)
        post_imag = np.ones((grid_num, grid_num))
        vmag=None
        imag_wavelength=None
        contrast_str=None
        imag_calc=None
        
    #######################################################################################
    ## Total
    #######################################################################################
    tot_list = hlp.tot_list(rv_list, astro_list, imag_list, num_points)
    
    if run_imag and imag_calc=='exact':
        post_tot = hlp.post_single(tot_list, a_inds, m_inds, grid_num)
        
    else:
        post_tot = hlp.post_tot_approx_imag(tot_list, post_imag, a_inds, m_inds, grid_num)
    #######################################################################################
    #######################################################################################
    
    ##
    end_time = time.time()
    ##
    if verbose:
        print('{:.0e} points ran for {} in {:.2f} seconds.'.format(num_points, star_name, end_time-start_time))
    
    if 'proc' in save:
        ls.save_processed(star_name, m_star, d_star,
                          run_rv, run_astro, run_imag, 
                          post_tot, post_rv, post_astro, post_imag,
                          prior, a_lim, m_lim, outdir=outdir, verbose=verbose)
    if 'raw' in save:
        if imag_calc=='approx':
            imag_data = post_imag
        else: # If imag_calc='exact' or anything else
            imag_data = imag_list
        
        ls.save_raw(star_name, m_star, d_star,
                    run_rv, run_astro, run_imag,
                    tot_list, rv_list, astro_list, imag_data,
                    prior, vmag, imag_wavelength, contrast_str,
                    a_list, m_list, a_inds, m_inds, a_lim, m_lim, 
                    imag_calc=imag_calc, outdir=outdir, 
                    verbose=verbose)     
    return
    
    
def plot(args):
    """
    Plot the content of loaded arrays

    Arguments:
        args: Command line arguments, including:
    
        config_path (str): Path to configuration file
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
    # Even though results are saved, config file still needed for a few parameters
    config_path = args.config
    cm = load_module_from_file(config_path)

    # Check if scatter_plot and outdir are provided in config. Otherwise set to defaults
    scatter_plot, outdir = set_values(config_path, ['scatter_plot', 'outdir'], [None, ''])

    star_name, m_star, d_star,\
    run_rv, run_astro, run_imag,\
    post_tot, post_rv, post_astro, post_imag,\
    prior, a_lim, m_lim = ls.load(args.read_file_path, cm.grid_num, args.verbose)

    if "2d" in args.type:
        plotter.joint_plot(star_name, m_star, d_star,
                           cm.run_rv, cm.run_astro, cm.run_imag,
                           post_tot, post_rv, post_astro, post_imag,
                           a_lim, m_lim,
                           scatter_plot=scatter_plot, period_lines=False,
                           outdir=outdir, verbose=args.verbose)

    if "1d" in args.type:
        plotter.plot_1d(star_name, post_tot, a_lim, m_lim, outdir=outdir)

    return
    
    
def lims(args):
    """
    Load in saved arrays and calculate 1D a and m bounds.
    Print out and return these bounds.
    
    Arguments:
        args: Command line arguments, including:

        config_path (str): Path to configuration file
        read_file_path (str): Path to saved outputs.
        grid_num (int): Shape of square posterior arrays.

    Returns:
        bounds (list of tuples): [(a1, a2),(m1, m2)] giving 95%
                                 confidence intervals
    """
    config_path = args.config
    cm = load_module_from_file(config_path)

    star_name, m_star, d_star,\
    run_rv, run_astro, run_imag,\
    post_tot, post_rv, post_astro, post_imag,\
    prior, a_lim, m_lim = ls.load(args.read_file_path, cm.grid_num, args.verbose)
    
    # bounds is the final answer: [range of 2σ a, range of 2σ m].
    # twosig_inds contains the indices corresponding to bounds. That is, where the CDF reaches the upper and lower values associated with the 95% confidence interval.
    bounds, twosig_inds = hlp.bounds_1D(post_tot, [m_lim, a_lim], 2)
        
    # Print out the 2-sigma boundaries (bounds) for the joint posterior
    # twosig_levels is a list of 2 floats: the 2sigma probs for a and m such that 95% of the prob is contained within the interval twosig_inds[i]
    print('a_lim = ', bounds[0], ' AU')
    print('m_lim = ', bounds[1], ' M_J')
    
    return bounds
    
def all(args):
    """
    Run the run, plot, and lims commands sequentially.
    """
    
    run(args)
    plot(args)
    lims(args)
    
    return
    
    
    
def load_module_from_file(config_path):
    """
    Adapted from radvel
    Loads a python module from the path of the corresponding file.
    Args:
        config_path (str): path to configuration file, 
                           e.g. "ethraid/config_files/default.py"
    Returns:
        A valid module object
    Raises:
        ImportError: when the module can't be loaded
        FileNotFoundError: when module_path doesn't exist
    """
    abs_path = os.path.abspath(config_path)
    
    if sys.version_info[0] == 3 and sys.version_info[1] >= 5:
        import importlib.util
        spec = importlib.util.spec_from_file_location(config_path, abs_path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        
    elif sys.version_info[0] == 3 and sys.version_info[1] < 5:
        import importlib.machinery
        loader = importlib.machinery.SourceFileLoader(config_path, abs_path)
        module = loader.load_module()
        
    elif sys.version_info[0] == 2:
        import imp
        module = imp.load_source(config_path, abs_path)

    return module
    
    
def set_values(config_path, param_names, default_values):
    """
    Check if each of a set of parameters is defined in the 
    configuration file. If it is, return the defined value. If
    it isn't, then return the provided default value.
    
    Arguments:
        config_path (str): Path to configuration file
        param_names (list of str): List of parameter names
        defaults_values (list): Default values of params in param_names.
                                Must have same length as param_names.
    
    Returns:
        param_values (tuple): Parameter values, whether default or given
                              in config file
    """
    
    if len(param_names) != len(default_values):
        raise Exception('Error: param_names and defaults_values must have the same length')
    
    if config_path == None: # If no config provided, then return default values
        param_values = default_values
    
    else:
    
        config_module = load_module_from_file(config_path)
        param_values = []
        for i in range(len(param_names)): # For each parameter name given
            param_name = param_names[i]
            default_value = default_values[i]
        
            try: # Try assigning the parameter value to a variable
                param_value = eval('config_module.{}'.format(param_name))
    
            except Exception as err: # If there is no such value, use the default instead
                param_value = default_value
            
            param_values.append(param_value)

        if len(param_values)==1:
            return param_values[0]
        
        else:
            param_values = tuple(param_values)
        
    return param_values
    
    

