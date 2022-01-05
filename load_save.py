import h5py
import numpy as np

import helper_functions_general as hlp


def load(read_file, grid_num, extension='posts/'):
    
    
    post_file_path = extension+read_file+'.h5'
    
    print('Reading posterior in from '+post_file_path)
    
    post_file = h5py.File(post_file_path, 'r')

    rv_list = np.array(post_file.get('rv_list')) # Probabilities associated with RV models
    astro_list = np.array(post_file.get('astro_list')) # Probabilities of astro models
    no_astro = np.array(post_file.get('no_astro')) # Bool indicating presence of astro data
    a_list = np.array(post_file.get('a_list')) # Semi-major axis values
    m_list = np.array(post_file.get('m_list')) # Companion mass values
    a_lim = np.array(post_file.get('a_lim'))
    m_lim = np.array(post_file.get('m_lim'))
    min_a, min_m = np.array(post_file.get('min_vals'))
    
    # Calculate indices using provided grid_num
    a_bins = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
    m_bins = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)

    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)
    
    min_index_m = int(np.ceil(hlp.value2index(min_m, (0, grid_num-1), m_lim)))
    min_index_a = int(np.ceil(hlp.value2index(min_a, (0, grid_num-1), a_lim)))
    
    
    prior_array = np.ones((grid_num, grid_num))
    prior_array[0:min_index_m, :] = 0
    prior_array[:, 0:min_index_a] = 0
    
    ##########################################
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
    
    return post_tot, post_rv, post_astro, a_lim, m_lim, min_a, min_m


def save(rv_list, astro_list, no_astro, a_list, m_list,
         a_lim, m_lim, min_a, min_m, write_file, extension='posts/'):
    
    post_file_path = extension+write_file+'.h5'
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
    
    return
    
    