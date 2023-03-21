import os
import numpy as np
import h5py

from ethraid.compiled import helper_functions_general as hlp
from ethraid import helper_functions_imaging as hlp_imag


def load(read_file_path, grid_num=None, verbose=False):
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
    
    ## Read in items that both raw and processed files have
    no_astro = np.array(post_file.get('no_astro')) # Bool indicating presence of astro data
    no_imag = np.array(post_file.get('no_imag')) # Bool indicating presence of astro data
    post_imag = np.array(post_file.get('post_imag')) # 2D array of probabilities from imaging
    

    a_lim = np.array(post_file.get('a_lim')) # Limits over which a is sampled
    m_lim = np.array(post_file.get('m_lim')) # Limits over which m is sampled
    
    data_type = post_file.get('data_type').asstr()[()]
    
    if data_type=='processed':
        
        if verbose:
            print('Searching for processed arrays')
            
            if 'grid_num' is not None:
                print("load_save.load: Argument 'grid_num' has no effect when loading processed arrays.\n"
                      "                Only raw arrays can be reshaped.")
        
        post_tot = np.array(post_file.get('post_tot')) # Marginalized probabity array associated with RV/astro models
        post_rv = np.array(post_file.get('post_rv')) # Marginalized probabity array associated with RV models
        post_astro = np.array(post_file.get('post_astro')) # Marginalized probabity array associated with astrometry models
        
        grid_num = np.shape(post_tot)[0]
        
        
    elif data_type=='raw':
        
        if verbose:
            print('Searching for raw arrays')

        assert grid_num is not None, "To load raw arrays, grid_num must be provided"
        
        rv_list = np.array(post_file.get('rv_list')) # Probability list associated with RV models
        astro_list = np.array(post_file.get('astro_list')) # Probability list of astro models
        a_list = np.array(post_file.get('a_list')) # Semi-major axis values
        m_list = np.array(post_file.get('m_list')) # Companion mass values
        
        try:
            vmag = np.array(post_file.get('vmag'))
            imag_wavelength = np.array(post_file.get('imag_wavelength'))
            contrast_str = post_file.get('contrast_str').asstr()[()]
            
        except Exception as err:
            vmag=None
            imag_wavelength=None
            contrast_str=None
            
            # Only print message if user provided imaging data. If they did not, then irrelevant
            if verbose and not no_imag:
                print("load_save.load: Error loading imaging parameters. \n " \
                      "               post_imag cannot be reshaped to a new grid_num.")
    
        # Calculate indices using provided grid_num
        a_bins = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
        m_bins = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)

        a_inds = np.digitize(a_list, bins = a_bins)
        m_inds = np.digitize(m_list, bins = m_bins)
    
        ##########################################
        # If no astrometry data, create a uniform array
        if no_astro:
            num_points = len(rv_list)
            astro_list = np.ones(num_points)
            post_astro = np.ones((grid_num, grid_num))
            
            if verbose:
                print("load_save.load: No astrometry data provided. \n" 
                      "                Bounds will be based on RVs only.")
        
        else:                                       
            post_astro = hlp.post_single(astro_list, a_inds, m_inds, grid_num)
        
        # post_imag is already shaped (see comments in helper_functions_imaging.py), so if you want to shape rv_list and astro_list to some other dimensions, post_imag must be recalculated to match that shape
        if np.shape(post_imag)[0] != grid_num:
            post_imag = hlp_imag.imag_array(d_star, vmag, imag_wavelength, 
                                            contrast_str, a_lim, m_lim, grid_num)
            
        post_rv = hlp.post_single(rv_list, a_inds, m_inds, grid_num)
        post_tot = hlp.post_tot(rv_list, astro_list, post_imag, grid_num, a_inds, m_inds)
    
    return star_name, m_star, d_star, post_tot, post_rv, post_astro, post_imag, grid_num, a_lim, m_lim


def save_raw(star_name, m_star, d_star, 
             rv_list, astro_list, post_imag, 
             no_astro, no_imag,
             vmag, imag_wavelength, contrast_str,
             a_list, m_list, a_lim, m_lim, 
             outdir='', verbose=False):
         
         """
         Saves raw 1D probability arrays to a specified h5py file.
         Note that you don't have to specify grid_num in order to save the arrays.
         You can save the raw arrays and later use load() to load them back in,
         specify grid_num then, and form them to whatever shape you want.
    
         Arguments:
             star_name (str): Name of star (does not need to be official)
             m_star (float, M_jup): Mass of host star
             d_star (float, AU): Distance from Earth to host star
             rv_list (list floats): RV data likelihoods conditioned
                                    on full orbit models
             astro_list (list of floats): Astrometry data likelihoods 
                                          conditioned on full orbit models
             post_imag (2D array of floats): Model probabilities given imaging 
                                             data only,marginalized over all 
                                             orbital parameters except a and m.
                                             NOTE that post_imag is a 2D array,
                                             not a list.
             no_astro (bool): Is astrometry data included?
             no_imag (bool): Is imaging data included?
             
             vmag (mag): Apparent V-band magnitude of host star
             imag_wavelength (μm): Wavelength of imaging observations
             contrast_str (dataframe or dict, 
                          columns of 'ang_sep' (arcseconds) 
                          and 'delta_mag' (mag)
                          ): Ordered pairs of angular separation and Δmag.

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
        
         save_dir = os.path.join(outdir, 'results/{}/'.format(star_name)) # Each star gets its own folder
         os.makedirs(save_dir, exist_ok=True)
         
         
         post_file_path = save_dir+star_name+'_raw.h5'
         post_file = h5py.File(post_file_path, 'w')
         
         # Save the un-binned values and arrays in case you want to use a different grid_num later
         post_file.create_dataset('star_name', data=star_name)
         post_file.create_dataset('m_star', data=m_star)
         post_file.create_dataset('d_star', data=d_star)
         
         post_file.create_dataset('rv_list', data=rv_list)
         post_file.create_dataset('astro_list', data=astro_list)
         post_file.create_dataset('post_imag', data=post_imag)
         post_file.create_dataset('no_astro', data=no_astro)
         post_file.create_dataset('no_imag', data=no_imag)
         
         try:
             # Save key imaging parameters only for the raw data, because if you want to reshape the raw 1d RV and astrometry arrays, you will have to recalculate the 2d imag_post array to match.
             post_file.create_dataset('vmag', data=vmag)
             post_file.create_dataset('imag_wavelength', data=imag_wavelength)
             post_file.create_dataset('contrast_str', data=contrast_str)
        
         except Exception as err:
             # Only print message if user provided imaging data. If they did not, then irrelevant
             if verbose and not no_imag:
                 print('load_save.save_raw: No imaging params saved. post_imag cannot be reshaped to a new grid_num.')
         
         post_file.create_dataset('a_list', data=a_list)
         post_file.create_dataset('m_list', data=m_list)
         
         post_file.create_dataset('a_lim', data=a_lim)
         post_file.create_dataset('m_lim', data=m_lim)
         
         post_file.create_dataset('data_type', data='raw')
         
         post_file.close()
         
         print('Posterior file saved to '+post_file_path)
         
         return
    
def save_processed(star_name, m_star, d_star, post_tot, 
                   post_rv, post_astro, post_imag,
                   no_astro, no_imag, a_lim, m_lim, outdir=''):
         
         """
         Saves shaped 2D probability arrays to a specified h5py file.
    
         Arguments:
             star_name (str): Name of star (does not need to be official)
             m_star (float, M_jup): Mass of host star
             d_star (float, AU): Distance from Earth to host star
             post_rv (2D array of floats): RV data likelihoods conditioned
                                           on full orbit models. Shaped into
                                           a 2D array, with model probabilities
                                           binned based on their (a,m) values
             post_astro (2D array of floats): Astrometry data likelihoods 
                                              conditioned on full orbit models. 
                                              Shaped into a 2D array, with 
                                              model probabilities binned based 
                                              on their (a,m) values
             post_imag (2D array of floats): Model probabilities given imaging 
                                             data only,marginalized over all 
                                             orbital parameters except a and m.
             no_astro (bool): Is astrometry data included?
             no_imag (bool): Is imaging data included?
                   
             a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                          in the form (a_min, a_max)
             m_lim (tuple of floats, M_jup): Mass limits as (m_min, m_max)
         
         Returns:
             None
         """
        
         save_dir = os.path.join(outdir, 'results/{}/'.format(star_name)) # Each star gets its own folder
         os.makedirs(save_dir, exist_ok=True)
         
         
         post_file_path = save_dir+star_name+'_processed.h5'
         post_file = h5py.File(post_file_path, 'w')
         
         # Save the un-binned values and arrays in case you want to use a different grid_num later
         post_file.create_dataset('star_name', data=star_name)
         post_file.create_dataset('m_star', data=m_star)
         post_file.create_dataset('d_star', data=d_star)
         
         post_file.create_dataset('post_tot', data=post_tot)
         post_file.create_dataset('post_rv', data=post_rv)
         post_file.create_dataset('post_astro', data=post_astro)
         post_file.create_dataset('post_imag', data=post_imag)
         post_file.create_dataset('no_astro', data=no_astro)
         post_file.create_dataset('no_imag', data=no_imag)
         
         post_file.create_dataset('a_lim', data=a_lim)
         post_file.create_dataset('m_lim', data=m_lim)
         
         post_file.create_dataset('data_type', data='processed')
         
         post_file.close()
         
         print('Posterior file saved to '+post_file_path)
         
         return