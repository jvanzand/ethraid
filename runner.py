import astropy.constants as c
import numpy as np
import time
import radvel as rv

import os
import sys

path = os.getcwd()
sys.path.append(path+'/trends')

import plotter
import system_params as sp

#########################
import helper_functions_general as hlp
import helper_functions_rv as hlp_rv
import helper_functions_astro as hlp_astro
import helper_functions_imaging as hlp_imag

import load_save as ls
#########################


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
M_earth = 5.972167867791379e+27

# General (necessary)
# RV
# Astrometry
# Imaging
# General (optional)
def run(star_name, m_star, d_star, 
        gammadot, gammadot_err, gammaddot, gammaddot_err, rv_baseline, rv_epoch, 
        delta_mu, delta_mu_err,
        vmag=None, imag_wavelength=None, contrast_str=None, 
        scatter_tuple=None, num_points=1e6, grid_num=100, save=True, plot=True, read_file_path=None):
        
    """
    Primary function to run trend code.
        
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
    
    Returns:
        None
    """

    # If no data to read in, calculate new arrays
    if read_file_path is None:
        
        min_per = rv_baseline*0.7
        # min_m = M_earth/M_jup
        min_m = 1
        m_star_Ms = m_star * M_jup/M_sun
        # Finally, the minimum semi-major axis is the one where period is smallest and companion mass is smallest too. If companion mass were larger at the same period, the companion would have to be farther away. Same for larger period at fixed mass.
        min_a = rv.utils.semi_major_axis(min_per, (m_star_Ms + min_m*(M_jup/M_sun)))

        print('Min sampling m is: ', min_m)
        print('Min sampling a is: ', min_a)

        max_a = 4e2
        max_m = 5e2
        
        # General
        a_lim = (min_a, max_a)
        m_lim = (min_m, max_m)

        num_points = int(num_points)
        # np.set_printoptions(threshold=np.inf)
    
        
        a_list, m_list, per_list, e_list, i_list,\
        om_list, M_anom_0_list, a_inds, m_inds = hlp.make_arrays(m_star, a_lim, m_lim,\
                                                                grid_num, num_points)

        print('made arrays')
        ##
        start_time = time.time()
        ##

        # Some targets aren't in the Hip/Gaia catalog, so we can't make the astrometry posterior for them.
        no_astro = False
        try:
            # Use a negative dmu value to run without astrometry
            if delta_mu < 0:
                raise ValueError()
            astro_list = hlp_astro.astro_list(a_list, m_list, e_list, i_list, 
                                              om_list, M_anom_0_list, per_list,
                                              m_star, d_star, delta_mu, delta_mu_err)                     
                                 
            post_astro = np.array(hlp.post_single(astro_list, a_inds, m_inds, grid_num))

        except Exception as err:
            astro_list = np.ones(num_points)
            post_astro = np.zeros((grid_num, grid_num))
            no_astro = True
            print('No astrometry data provided. Bounds will be based on RVs only.')
    


        rv_list = hlp_rv.rv_list(a_list, m_list, e_list, i_list, om_list, M_anom_0_list,
                                per_list, m_star, rv_epoch,
                                gammadot, gammadot_err, gammaddot, gammaddot_err)
                                
        post_rv = hlp.post_single(rv_list, a_inds, m_inds, grid_num)
        
        post_imag = hlp_imag.imag_array(d_star, vmag, imag_wavelength, contrast_str, a_lim, m_lim, grid_num)
        
        post_tot = hlp.post_tot(rv_list, astro_list, post_imag, grid_num, a_inds, m_inds)
        
        ##
        end_time = time.time()
        ##
        print('{:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time-start_time))

    
        if save==True:
            ls.save(star_name, rv_list, astro_list, no_astro, post_imag, a_list, m_list,
                    a_lim, m_lim, min_a, min_m)
    
    # Otherwise, load in existing data:
    else:
        post_tot, post_rv, post_astro, post_imag, a_lim, m_lim, min_a, min_m = ls.load(read_file_path, grid_num)
        
        
    if plot==True:
        plotter.joint_plot(star_name, m_star, post_tot, post_rv, post_astro, post_imag, grid_num, 
                a_lim, m_lim, scatter_tuple=scatter_tuple, period_lines = False)
    
    return
    
    
    

if __name__ == "__main__":
    
    # import pandas as pd
    # contrast_curve = pd.read_csv('data/EDG_raw_curves/191939_832_sensitivity.dat',
    #                               skiprows=29,
    #                               delimiter=' ',
    #                               header=None)
    # new_header = ['ang_sep', 'delta_mag']
    # contrast_curve.columns = new_header
    # contrast_curve.to_csv('data/EDG_clean_curves/191939_832_clean.csv', index=False)
    
    run(*sp.params_191939_old, num_points=1e6, grid_num=100, plot=True, read_file_path=None)
    # 'results/post_arrays/T001174.h5')
    #'results/post_arrays/12572.h5')
    # run(*sp.params_synth, num_points=1e6, grid_num=100, save=False, plot=True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

