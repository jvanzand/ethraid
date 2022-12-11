import os
import numpy as np
import h5py

from trends import helper_functions_general as hlp


def load(read_file_path, grid_num):
    
    
    print('Reading posterior in from '+read_file_path)
    
    post_file = h5py.File(read_file_path, 'r')
    
    star_name = post_file.get('star_name').asstr()[()] # Special treatment of H5py strings    
    m_star = np.array(post_file.get('m_star'))
    d_star = np.array(post_file.get('d_star'))

    rv_list = np.array(post_file.get('rv_list')) # Probabilities associated with RV models
    astro_list = np.array(post_file.get('astro_list')) # Probabilities of astro models
    no_astro = np.array(post_file.get('no_astro')) # Bool indicating presence of astro data
    post_imag = np.array(post_file.get('post_imag')) # 2D array of probabilities from imaging
    # post_imag = np.ones((grid_num, grid_num))
    a_list = np.array(post_file.get('a_list')) # Semi-major axis values
    m_list = np.array(post_file.get('m_list')) # Companion mass values
    a_lim = np.array(post_file.get('a_lim')) # Limits over which a is sampled
    m_lim = np.array(post_file.get('m_lim')) # Limits over which m is sampled
    min_a, min_m = np.array(post_file.get('min_vals'))
    
    # Calculate indices using provided grid_num
    a_bins = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
    m_bins = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)

    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)
    
    ##########################################
    if no_astro:
        num_points = len(rv_list)
        astro_list = np.ones(num_points)
        post_astro = np.zeros((grid_num, grid_num))
        
        grid_pad = int(np.round(grid_num/15))
        post_astro = np.pad(post_astro, [(grid_pad, 0), (grid_pad, 0)])
        print('No astrometry data provided. Bounds will be based on RVs only.')
        
    else:                                       
        post_astro = hlp.post_single(astro_list, a_inds, m_inds, grid_num)
        # post_astro = post_astro/post_astro.sum()
    
    post_rv = hlp.post_single(rv_list, a_inds, m_inds, grid_num)
    # post_rv = post_rv/post_rv.sum()

    post_tot = hlp.post_tot(rv_list, astro_list, post_imag, grid_num, a_inds, m_inds)
    # post_tot = post_tot/post_tot.sum()
    
    return star_name, m_star, d_star, post_tot, post_rv, post_astro, post_imag, a_lim, m_lim, min_a, min_m


def save(star_name, m_star, d_star, rv_list, astro_list, no_astro, post_imag,
         a_list, m_list, a_lim, m_lim, min_a, min_m):
        
    save_dir = 'results/post_arrays/' # Arrays for all stars go in one folder
    os.makedirs(save_dir, exist_ok=True)
    
    
    post_file_path = save_dir+star_name+'.h5'
    post_file = h5py.File(post_file_path, 'w')
    
    # Save the un-binned values and arrays in case you want to use a different grid_num later
    post_file.create_dataset('star_name', data=star_name)
    post_file.create_dataset('m_star', data=m_star)
    post_file.create_dataset('d_star', data=d_star)
    
    post_file.create_dataset('rv_list', data=rv_list)
    post_file.create_dataset('astro_list', data=astro_list)
    post_file.create_dataset('no_astro', data=no_astro)
    post_file.create_dataset('post_imag', data=post_imag)
    
    post_file.create_dataset('a_list', data=a_list)
    post_file.create_dataset('m_list', data=m_list)
    
    post_file.create_dataset('a_lim', data=a_lim)
    post_file.create_dataset('m_lim', data=m_lim)
    
    post_file.create_dataset('min_vals', data=(min_a, min_m))
    
    post_file.close()
    print('Posterior file saved to '+post_file_path)
    
    return
    
    