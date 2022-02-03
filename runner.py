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

import load_save as ls
#########################


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
M_earth = 5.972167867791379e+27

# params_star = (m_star, distance(cm), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), rv_range, rv_epoch, delta_mu, delta_mu_err)
def run(star_name, m_star, d_star, gammadot, gammadot_err, gammaddot, gammaddot_err,
        rv_baseline, rv_epoch, delta_mu, delta_mu_err,
        num_points=1e6, grid_num=100, save=True, plot=True, 
        read_file_path=None):
        
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
    """

    # If no data to read in, calculate new arrays
    if read_file_path is None:
        
        ##### Determination of minimum mass constraint. Kind of a headache but could be useful.
        ##### Replacing for now with blanket m_min = 1 Earth mass
        # # RV 0 point is arbitrary, so the min. K amp. is half of the 'peak to trough' of the RVs.
        # # min_K = rv_range/2 # This is prone to errors because it relies on the individual end points.
        #
        # m_star_Ms = m_star * M_jup/M_sun
        # min_per = rv_baseline
        # # This still has the issue that the end of the timeseries might not be the most extreme point, namely if there is enough curvature to reach a maximum and then come back down (or up). Fortunately, this just means I would be sampling masses that were too small, rather than cutting any out. Still could use some refinement.
        # min_K = abs(0.5*(gammadot*rv_baseline + 0.5*gammaddot*rv_baseline**2))
        #
        #
        # # The argument for min_m is this: suppose period is min_per, which is where m can be smallest. What is the smallest that m can be? It can be so small that the current max RV is the highest the RV will ever be. Then m is so small that 1) it attains the lowest possible K and 2) it does so as close in as possible (any farther out would mean a larger planet). NOTE: this happens at e=0, which is NOT concordant with minimum m. See below for idea on how to fix this.
        # # One way to make this more general would be to choose an eccentricity that is ~2σ from the most likely value. Roughly speaking, we could then say we were only cutting out 2σ discrepant models for min_m.
        # min_m = rv.utils.Msini(min_K, min_per, m_star_Ms, e=0, Msini_units='jupiter')
        ####################################################################################
        
        min_per = rv_baseline
        # min_m = M_earth/M_jup
        min_m = 1
        m_star_Ms = m_star * M_jup/M_sun
        # Finally, the minimum semi-major axis is the one where period is smallest and companion mass is smallest too. If companion mass were larger at the same period, the companion would have to be farther away. Same for larger period at fixed mass.
        min_a = rv.utils.semi_major_axis(min_per, (m_star_Ms + min_m*(M_jup/M_sun)))

        ###################################
        # # Experimental: revised min_m using gdot. The idea is to find the smallest mass that could produce the observed gdot at the known minimum period. This is not completely right because it uses min_K to get min_m, min_m to get min_a, and then min_a to get a new value for min_m.
        #
        # min_m = (gammadot*100/(24*3600))*((min_per*24*3600)/6.283185)**2*(m_star*M_sun)/(min_a*14959787070000.0) / M_jup
        ###################################

        print('Min sampling m is: ', min_m)
        print('Min sampling a is: ', min_a)

        # Sampling limits for a and m.
        # # 191939
        # # min_a = 0.5
        # # min_m = 0.5
        # a_lim = (0.8*min_a, 5e1)
        # m_lim = (0.8*min_m, 1e2)
    
        max_a = 1e2
        max_m = 1e3
        
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
                raise ValueError('delta_mu is less than 0.')
            astro_list = hlp_astro.astro_list(a_list, m_list, e_list, i_list, 
                                              om_list, M_anom_0_list, per_list,
                                              m_star, d_star, delta_mu, delta_mu_err)                     
                                 
            post_astro = np.array(hlp.prob_array(astro_list, a_inds, m_inds, grid_num))
            post_astro = post_astro/post_astro.sum()

        except Exception as err:
            print(err)
            astro_list = np.ones(num_points)
            post_astro = np.ones((grid_num, grid_num))
            no_astro = True
            print('No astrometry data provided. Bounds will be based on RVs only.')
    


        rv_list = hlp_rv.rv_list(a_list, m_list, e_list, i_list, om_list, M_anom_0_list,
                                per_list, m_star, rv_epoch,
                                gammadot, gammadot_err, gammaddot, gammaddot_err)
                                
        post_rv = np.array(hlp.prob_array(rv_list, a_inds, m_inds, grid_num))
        post_rv = post_rv/post_rv.sum()
        
        post_tot = np.array(hlp.post_tot(rv_list, astro_list, grid_num, a_inds, m_inds))
        post_tot = post_tot/post_tot.sum()

        ##
        end_time = time.time()
        ##
        print('{:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time-start_time))

    
        if save==True:
            ls.save(star_name, rv_list, astro_list, no_astro, a_list, m_list,
                    a_lim, m_lim, min_a, min_m)
    
    # Otherwise, load in existing data:
    else:
        post_tot, post_rv, post_astro, a_lim, m_lim, min_a, min_m = ls.load(read_file_path, grid_num)
        
        
    if plot==True:
        plotter.joint_plot(star_name, m_star, post_tot, post_rv, post_astro, grid_num, 
                a_lim, m_lim, period_lines = False)
    
    return
    
    
    

if __name__ == "__main__":
    
    run(*sp.params_hd182488, num_points=1e7, grid_num=70, plot=True, read_file_path='results/post_arrays/hd182488.h5')
    #'results/post_arrays/12572.h5')
    # run(*sp.params_synth, num_points=1e6, grid_num=100, save=False, plot=True)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

