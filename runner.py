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

# params_star = (m_star, distance(cm), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), max_rv of residuals, rv_epoch, delta_mu, delta_mu_err)
def run(m_star, d_star, gammadot, gammadot_err, gammaddot, gammaddot_err,
        rv_baseline, max_rv, rv_epoch, delta_mu, delta_mu_err,
        num_points=1e6, grid_num=100, save=True, plot=True, 
        read_file=None, write_file=None):
    
    # m_star, d_star, gammadot, gammadot_err, gammaddot, gammaddot_err,\
    #         rv_baseline, max_rv, rv_epoch, delta_mu, delta_mu_err = params


    # min_per is 4xbaseline for 191939 because we see ~no curvature yet.
    # min_per = 4*rv_baseline
    min_per = rv_baseline
    min_K = max_rv
    
    m_star_Ms = m_star * M_jup/M_sun
    # One way to make this more general would be to choose an eccentricity that is ~2σ from the
    # Most likely value. Roughly speaking, we could then say we were only cutting out 2σ discrepant models.
    min_m = rv.utils.Msini(min_K, min_per, m_star_Ms, e=0, Msini_units='jupiter')
    min_a = rv.utils.semi_major_axis(min_per, (m_star_Ms + min_m*(M_jup/M_sun)))

    ###################################
    # # Experimental: revised min_m using gdot. The idea is to find the smallest mass that could produce the observed gdot at the known minimum period. This is not completely right because it uses min_K to get min_m, min_m to get min_a, and then min_a to get a new value for min_m.
    #
    # min_m = (gammadot*100/(24*3600))*((min_per*24*3600)/6.283185)**2*(m_star*M_sun)/(min_a*14959787070000.0) / M_jup
    ###################################

    print('Min m is: ', min_m)
    print('Min a is: ', min_a)

    # Sampling limits for a and m. Note that if the min_a or min_m parameters fall outside these bounds, the plot will look weird. I can modify later to throw an error, but it's mostly visual.
    # # 191939
    # # min_a = 0.5
    # # min_m = 0.5
    # a_lim = (0.8*min_a, 5e1)
    # m_lim = (0.8*min_m, 1e2)
    
    # General
    a_lim = (0.8*min_a, 1e2)
    m_lim = (0.8*min_m, 2e2)
    print(a_lim[0], min_a)

    num_points = int(num_points)
    np.set_printoptions(threshold=np.inf)
    
    
    # If you are loading in existing data:
    if read_file is not None:
        
        a_inds, m_inds, min_index_m, min_index_a, prior_array = ls.load(read_file, 
                                                                        extension='posts/')
        
        if no_astro:
            num_points = len(rv_list)
            astro_list = np.ones(num_points)
            post_astro = np.ones((grid_num, grid_num))
            print('No astrometry data provided. Bounds will be based on RVs only.')
            
        else:                                       
            post_astro = hlp.prob_array(astro_list, a_inds, m_inds, grid_num) * prior_array
            post_astro = post_astro/post_astro.sum()
        
        post_rv = hlp.prob_array(rv_list, a_inds, m_inds, grid_num) * prior_array
        post_rv = post_rv/post_rv.sum()

        post_tot = hlp.post_tot(rv_list, astro_list, grid_num, a_inds, m_inds) * prior_array
        post_tot = post_tot/post_tot.sum()
        
    # Otherwise, calculate new arrays
    else:
        
        a_list, m_list, per_list, e_list, i_list,\
        om_list, M_anom_0_list, a_inds, m_inds = hlp.make_arrays(m_star, a_lim, m_lim, rv_epoch,\
                                                                grid_num, num_points)

        print('made arrays')
        ##
        start_time = time.time()
        ##

        # Create an array with 1s in allowed regions and 0s in disallowed regions
        min_index_m = int(np.ceil(hlp.value2index(min_m, (0, grid_num-1), m_lim)))
        min_index_a = int(np.ceil(hlp.value2index(min_a, (0, grid_num-1), a_lim)))

        prior_array = np.ones((grid_num, grid_num))
        prior_array[0:min_index_m, :] = 0
        prior_array[:, 0:min_index_a] = 0

        # Some targets aren't in the Hip/Gaia catalog, so we can't make the astrometry posterior for them.
        no_astro = False
        try:
            astro_list = hlp_astro.astro_list(a_list, m_list, e_list, i_list, 
                                              om_list, M_anom_0_list, per_list,
                                              m_star, d_star, delta_mu, delta_mu_err)                     
                                 
            post_astro = hlp.prob_array(astro_list, a_inds, m_inds, grid_num) * prior_array
            post_astro = post_astro/post_astro.sum()

        except Exception as e:
            print(e)
            astro_list = np.ones(num_points)
            post_astro = np.ones((grid_num, grid_num))
            no_astro = True
            print('No astrometry data provided. Bounds will be based on RVs only.')
    


        rv_list = hlp_rv.rv_list(a_list, m_list, e_list, i_list, om_list, M_anom_0_list,
                                per_list, m_star, rv_epoch,
                                gammadot, gammadot_err, gammaddot, gammaddot_err)
                                
        post_rv = hlp.prob_array(rv_list, a_inds, m_inds, grid_num) * prior_array
        post_rv = post_rv/post_rv.sum()
        
        post_tot = hlp.post_tot(rv_list, astro_list, grid_num, a_inds, m_inds) * prior_array
        post_tot = post_tot/post_tot.sum()

        ##
        end_time = time.time()
        ##
        print('{:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time-start_time))

    
        if save==True:
            ls.save(rv_list, astro_list, no_astro, a_list, m_list,
                    a_lim, m_lim, min_a, min_m, write_file, extension='posts/')
        if plot==True:
            plotter.joint_plot(m_star, post_tot, post_rv, post_astro, grid_num, a_lim, m_lim, (min_a, min_m),
                    save_name='base', period_lines = False)
    
    return


if __name__ == "__main__":
    
    run(*sp.params_191939, num_points=1e6, grid_num=100, save=True, plot=True, read_file=None, write_file='base')

