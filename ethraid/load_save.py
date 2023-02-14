import os
import numpy as np
import h5py

from ethraid.compiled import helper_functions_general as hlp


def load(read_file_path, grid_num):
    """
    Loads probability arrays from a specified h5py file.
    
    Arguments:
        read_file_path (str): Path to saved data
        grid_num (int): Desired array shape
    
    Returns:
        star_name (str): Name of star (does not need to be official)
        m_star (float, M_jup): Mass of host star
        d_star (float, AU): Distance from Earth to host star
        post_tot (array of floats): Total posterior array, shape=(grid_num,grid_num)
        post_rv (array of floats): Model probabilities given RV data only, marginalized
                                   over all orbital parameters except a and m
        post_astro (array of floats): Model probabilities given astrometry data only,
                                      marginalized over all orbital parameters except 
                                      a and m
        post_imag (array of floats): Model probabilities given imaging data only, 
                                     marginalized over all orbital parameters 
                                     except a and m
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max)
        m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
    """
    
    
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
    
    # Calculate indices using provided grid_num
    a_bins = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
    m_bins = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)

    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)
    
    ##########################################
    if no_astro:
        num_points = len(rv_list)
        astro_list = np.ones(num_points)
        post_astro = np.ones((grid_num, grid_num))

        print('No astrometry data provided. Bounds will be based on RVs only.')
        
    else:                                       
        post_astro = hlp.post_single(astro_list, a_inds, m_inds, grid_num)
    post_rv = hlp.post_single(rv_list, a_inds, m_inds, grid_num)
    post_tot = hlp.post_tot(rv_list, astro_list, post_imag, grid_num, a_inds, m_inds)
    
    return star_name, m_star, d_star, post_tot, post_rv, post_astro, post_imag, a_lim, m_lim


def save(star_name, m_star, d_star, rv_list, astro_list, no_astro, post_imag,
         a_list, m_list, a_lim, m_lim, outdir=''):
         
         """
         Saves calculated probability arrays to a specified h5py file.
         Note that you don't have to specify grid_num in order to save the arrays.
         You can save the raw arrays and later use load() to load them back in,
         specify grid_num then, and form them to whatever shape you want.
    
         Arguments:
             read_file_path (str): Path to saved data
             grid_num (int): Desired array shape
    
         Arguments:
             star_name (str): Name of star (does not need to be official)
             m_star (float, M_jup): Mass of host star
             d_star (float, AU): Distance from Earth to host star
             rv_list (list floats): RV data likelihoods conditioned
                                    on full orbit models
             astro_list (list of floats): Astrometry data likelihoods 
                                          conditioned on full orbit models
             no_astro (bool): Is astrometry data included?
             post_imag (array of floats): Model probabilities given imaging 
                                          data only,marginalized over all 
                                          orbital parameters except a and m.
                                          NOTE that post_imag is a 2D array,
                                          not a list.
             a_list (list of floats): Semi-major axes corresponding to above
                                      model probabilities
             m_list (list of floats): Companion masses corresponding to above
                                      model probabilities
             a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                          in the form (a_min, a_max)
             m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
         
         Returns:
             None
         """
        
         save_dir = outdir+'results/post_arrays/'# Arrays for all stars go in one folder
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
         
         post_file.close()
         print('Posterior file saved to '+post_file_path)
         
         return
    
    