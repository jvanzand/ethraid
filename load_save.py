import h5py
import numpy as np


def load(read_file, extension='posts/'):
    
    
    post_file_path = extension+read_file+'.h5'
    
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
    
    min_index_m = int(np.ceil(hlp.value2index(min_m, (0, grid_num-1), m_lim)))
    min_index_a = int(np.ceil(hlp.value2index(min_a, (0, grid_num-1), a_lim)))
    
    
    prior_array = np.ones((grid_num, grid_num))
    prior_array[0:min_index_m, :] = 0
    prior_array[:, 0:min_index_a] = 0
    
    return a_inds, m_inds, min_index_m, min_index_a, prior_array


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
    
    