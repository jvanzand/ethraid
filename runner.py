import matplotlib.pyplot as plt
import astropy.constants as c
import numpy as np
# from astropy.time import Time
from scipy.stats import loguniform, beta
import time
import h5py

import radvel as rv


from trends import helper_functions_wrapper as hlpw
import plotter
import system_params as sp


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
      

def run(read_file=None, write_file=None, num_points=1e6, grid_num=100, save=True):
    
    # rv_epoch is the epoch where DATA values of g_dot and g_ddot are computed. Taken from radvel setup file.
    m_star, d_star, gammadot, gammadot_err, gammaddot, gammaddot_err,\
            rv_baseline, max_rv, rv_epoch, delta_mu, delta_mu_err = sp.params_synth


    # min_per is 4xbaseline because we see ~no curvature yet.
    min_per = 4*rv_baseline
    # min_per = rv_baseline
    min_K = max_rv

    min_m = rv.utils.Msini(min_K, min_per, m_star, e=0, Msini_units='jupiter')
    min_a = rv.utils.semi_major_axis(min_per, (m_star + min_m*(M_jup/M_sun)))

    # # Experimental: revised min_m using gdot. The idea is to find the smallest mass that could produce the observed gdot at the known minimum period. This is not completely right because it uses min_K to get min_m, min_m to get min_a, and then min_a to get a new value for min_m.
    #
    # min_m = (gammadot*100/(24*3600))*((min_per*24*3600)/6.283185)**2*(m_star*M_sun)/(min_a*14959787070000.0) / M_jup

    print('Min m is: ', min_m)
    print('Min a is: ', min_a)

    # Sampling limits for a and m. Note that if the min_a or min_m parameters fall outside these bounds, the plot will look weird. I can modify later to throw an error, but it's mostly visual.
    # # 191939
    # # min_a = 0.5
    # # min_m = 0.5
    # a_lim = (0.8*min_a, 5e1)
    # m_lim = (0.8*min_m, 1e2)
    # # HD238894 - Paul Dalba's target
    # a_lim = (0.8*min_a, 32.1)
    # m_lim = (0.8*min_m, 3e2)
    # # HD6101
    # a_lim = (0.9*min_a, 6e1)
    # m_lim = (0.9*min_m, 1e5)
    # # HD91204
    # a_lim = (0.8*min_a, 8e1)
    # m_lim = (0.8*min_m, 1e5)
    # synthetic
    a_lim = (0.8*min_a, 5e1)
    m_lim = (0.8*min_m, 1e3)
    # # GL758
    # a_lim = (0.5*min_a, 2e2)
    # m_lim = (0.5*min_m, 4e2)
    # # HIP67246
    # a_lim = (0.5*min_a, 2e2)
    # m_lim = (0.5*min_m, 4e2)
    # # HIP63510
    # a_lim = (0.5*min_a, 2e2)
    # m_lim = (0.5*min_m, 1e3)
    # # 12572
    # a_lim = (0.5*min_a, 5e1)
    # m_lim = (0.5*min_m, 1e2)
    print(a_lim[0], min_a)

    num_points = int(num_points)
    np.set_printoptions(threshold=np.inf)
    
    
    # If you are loading in existing data
    if read_file is not None:
        
        # Load data
        post_file_path = 'posts/'+read_file+'.h5'
        print('Reading posterior in from '+post_file_path)
        post_file = h5py.File(post_file_path, 'r')

        rv_list = np.array(post_file.get('rv_list'))
        astro_list = np.array(post_file.get('astro_list'))
        no_astro = np.array(post_file.get('no_astro'))
        a_list = np.array(post_file.get('a_list'))
        m_list = np.array(post_file.get('m_list'))
        a_lim = np.array(post_file.get('a_lim'))
        m_lim = np.array(post_file.get('m_lim'))
        min_a, min_m = np.array(post_file.get('min_vals'))
        
        # Calculate indices using provided grid_num
        a_bins = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
        m_bins = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)

        a_inds = np.digitize(a_list, bins = a_bins)
        m_inds = np.digitize(m_list, bins = m_bins)
        
        min_index_m = int(np.ceil(hlpw.value2index(min_m, (0, grid_num-1), m_lim)))
        min_index_a = int(np.ceil(hlpw.value2index(min_a, (0, grid_num-1), a_lim)))
        
        
        prior_array = np.ones((grid_num, grid_num))
        prior_array[0:min_index_m, :] = 0
        prior_array[:, 0:min_index_a] = 0
        
        if no_astro:
            num_points = len(rv_list)
            astro_list = np.ones(num_points)
            post_astro = np.ones((grid_num, grid_num))
            print('No astrometry data provided. Bounds will be based on RVs only.')
            
        else:                                       
            post_astro = hlpw.prob_array(astro_list, a_inds, m_inds, grid_num) * prior_array
            post_astro = post_astro/post_astro.sum()
        
        post_rv = hlpw.prob_array(rv_list, a_inds, m_inds, grid_num) * prior_array
        post_rv = post_rv/post_rv.sum()

        post_tot = hlpw.post_tot(rv_list, astro_list, grid_num, a_inds, m_inds) * prior_array
        post_tot = post_tot/post_tot.sum()
        
    # Otherwise, calculate new arrays
    else:
        
        a_list, m_list, per_list, e_list, i_list, om_list, M_anom_0, E_anom_rv, a_inds, m_inds = \
                                                    hlpw.make_arrays(m_star, a_lim, m_lim, rv_epoch, grid_num, num_points)


        print('made arrays')

        ##
        start_time = time.time()
        ##

        # Create an array with 1s in allowed regions and 0s in disallowed regions
        min_index_m = int(np.ceil(hlpw.value2index(min_m, (0, grid_num-1), m_lim)))
        min_index_a = int(np.ceil(hlpw.value2index(min_a, (0, grid_num-1), a_lim)))

        prior_array = np.ones((grid_num, grid_num))
        prior_array[0:min_index_m, :] = 0
        prior_array[:, 0:min_index_a] = 0

        # plt.imshow(prior_array, origin='lower')
        # plt.show()

        # Some targets aren't in the Hip/Gaia catalog, so we can't make the astrometry posterior for them.
        no_astro = False
        try:
            astro_list = hlpw.astro_post(delta_mu, delta_mu_err, m_star, d_star, a_list,
                                         m_list, per_list, e_list, i_list, om_list,
                                         M_anom_0, num_points, grid_num)                      
                                 
            post_astro = hlpw.prob_array(astro_list, a_inds, m_inds, grid_num) * prior_array
            post_astro = post_astro/post_astro.sum()

        except:
            astro_list = np.ones(num_points)
            post_astro = np.ones((grid_num, grid_num))
            no_astro = True
            print('No astrometry data provided. Bounds will be based on RVs only.')
    


        rv_list = hlpw.rv_post(gammadot, gammadot_err, gammaddot, gammaddot_err, m_star, 
                                a_list, m_list, per_list, e_list, i_list, om_list, E_anom_rv, 
                                num_points, grid_num)
                                
        post_rv = hlpw.prob_array(rv_list, a_inds, m_inds, grid_num) * prior_array
        post_rv = post_rv/post_rv.sum()
        
        post_tot = hlpw.post_tot(rv_list, astro_list, grid_num, a_inds, m_inds) * prior_array
        post_tot = post_tot/post_tot.sum()

        ##
        end_time = time.time()
        ##
        print('{:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time-start_time))
    
        if save==True:
            post_file_path = 'posts/'+write_file+'.h5'
            post_file = h5py.File(post_file_path, 'w')
            
            # Save the un-binned arrays in case you want to use a different grid_num later
            post_file.create_dataset('rv_list', data=rv_list)
            post_file.create_dataset('astro_list', data=astro_list)
            post_file.create_dataset('no_astro', data=no_astro)
            
            post_file.create_dataset('a_list', data=a_list)
            post_file.create_dataset('m_list', data=m_list)
            
            post_file.create_dataset('a_lim', data=a_lim)
            post_file.create_dataset('m_lim', data=m_lim)
            
            post_file.create_dataset('min_vals', data=(min_a, min_m))
            
            post_file.close()
            print('Posterior file saved to '+post_file_path)
    
    return m_star, post_tot, post_rv, post_astro, grid_num, a_lim, m_lim, (min_a, min_m)

if __name__ == "__main__":
    
    
    m_star, post_tot, post_rv, post_astro, grid_num, a_lim, m_lim, (min_a, min_m) = \
            run(read_file=None, save=True, write_file='test_post', num_points=1e6, grid_num=100)


    plotter.joint_plot(m_star, post_tot, post_rv, post_astro, grid_num, a_lim, m_lim, (min_a, min_m), 
                        save_name='test', period_lines = False)




# plt.imsave('post_rv.png', post_rv, origin='lower')
# plt.imsave('post_astro.png', post_astro, origin='lower')
# plt.imsave('post_tot.png', post_tot, origin='lower')
#
# plt.imshow(post_rv, origin='lower')
# plt.show()
# plt.imshow(post_astro, origin='lower', cmap='jet')
# plt.show()
# plt.imshow(post_tot, origin='lower')
# plt.show()