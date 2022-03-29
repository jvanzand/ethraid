import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

import helper_functions_general as hlp


mamajek_table = pd.read_csv('data/mamajek.csv').rename(columns={'K_s':'K'})
baraffe_table = pd.read_csv('data/baraffe_table_4.csv')
bands= pd.read_csv('data/bands.csv').replace('K_s', 'K')
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value


def imag_array(d_star, vmag, imag_wavelength, contrast_str, a_lim, m_lim, grid_num):
    """
    Convert a contrast curve in (angular_separation, Δmag) space
    into (separation, mass) space. This requires finding a
    correspondence between brightness and mass, which I do with
    either the Mamajek table or the Baraffe (03) paper, depending
    on the mass regime.
    
    Arguments:
        d_star (float, pc): Distance to the host star
        vmag (float, mag): Apparent V-band magnitude of host star
        imag_wavelength (float, μm): Wavelength of imaging data in contrast_curve
        contrast_curve (str): Path to contrast curve with columns of 
                              'ang_sep' (arcseconds) and 'delta_mag' (mag)
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max).
        m_lim (tuple of floats, M_jup): Mass limits, (m_min, m_max).
        grid_num (int): Dimensions of (a,m) grid.
    
    returns:
        np_imag_array (2D array, dim=grid_numXgrid_num): Probability
                        array from imaging data. In this model,
                        imaging probabilities are either 1 or 0.
    """
    if vmag is None or imag_wavelength is None or contrast_str is None:
        return np.ones((grid_num, grid_num))
    
    # Get band of observations and host star abs mag in that band
    band_name, host_abs_Xmag = abs_Xmag(d_star, vmag, imag_wavelength)
    
    contrast_curve = pd.read_csv(contrast_str)
    
    if not set(['ang_sep', 'delta_mag']).issubset(contrast_curve.columns):
        raise Exception("The dataframe must contain columns 'ang_sep' and 'delta_mag'")
    # Objective is sep (AU) vs. mass (M_jup)
    ############## 1: Get sep. Make sure to convert d_star from au to pc
    seps = contrast_curve['ang_sep']*(d_star/pc_in_au)
    ##############
    
    ############## 2: Convert absolute mag to mass
    Xband_abs_mags = contrast_curve['delta_mag']+host_abs_Xmag
    
    # Start with Mamajek. We could cut super bright targets from the interpolation, but no need
    interp_df_mamajek = mamajek_table[['M_jup', band_name]]
    
    # It may be that Mamajek covers the whole contrast curve and we don't need Baraffe, but maybe not.
    # Find the dimmest mag in Mamajek and start selecting Baraffe mags after that. 
    max_mag = interp_df_mamajek[band_name].max()

    # Just like we took super bright entries from Mamajek, take even the dimmest from Baraffe to be safe.
    interp_df_baraffe = baraffe_table.query("{}>{}".format(band_name, max_mag))[['M_jup', band_name]]

    # Concatenate the two dfs above
    interp_df = pd.concat([interp_df_mamajek, interp_df_baraffe]).sort_values(by=band_name)
    
    
    companion_masses = mag_to_mass(interp_df[band_name], interp_df['M_jup'], Xband_abs_mags)
    ##############
    
    # Discrete correspondence between separation and mass. We get a continuous correspondence below
    a_m_contrast = pd.DataFrame({'M_jup':companion_masses, 'sep':seps}).reset_index(drop=True)
    
    # We might want to go out to sma values beyond what's given in the contrast curve.
    # In that case, conservatively estimate that the curve becomes flat after the last sma value
    # (It would actually continue to drop to lower masses, but increasingly slowly, so this is a good approximation)
    last_a_ind = a_m_contrast[a_m_contrast['sep'] == a_m_contrast['sep'].max()].index # Find largest a
    last_m = a_m_contrast.iloc[last_a_ind]['M_jup'] # Find m corresponding to largest a
    
    a_m_interp_fn = interp1d(a_m_contrast['sep'], a_m_contrast['M_jup'], 
                             bounds_error=False, fill_value=last_m)
    
    # import matplotlib.pyplot as plt
    # a_list = np.linspace(8, 100, 40)
    # plt.plot(a_list, a_m_interp_fn(a_list))
    # plt.show()
    # fdf
                             
    imag_array = np.ones((grid_num, grid_num))
    
    
    # Note that this is a coarse approximation. The "mass" over the whole ith grid square is just the conversion of the ith index to mass (and same for sma).
    for i in range(grid_num): # on mass axis
        m = hlp.index2value(i, (0, grid_num), m_lim)
        for j in range(grid_num): # on a axis
            a = hlp.index2value(j, (0, grid_num), a_lim)
            
            max_m = a_m_interp_fn(a)
            
            if m > max_m:
                imag_array[i, j] = 0
                
    np_imag_array = np.array(imag_array) / np.array(imag_array).sum()
    
    return np_imag_array


# An important question I might ask later: why do I jump straight to “imag_array” (100X100) rather than first creating “imag_list” (len=1e8) and forming it into an array, as I do for RVs and astrometry? Answer: Basically, it takes too long and it's not needed. The constraints imposed on imaging are fairly rough (in my model): at a given separation, rule out all masses greater than X. Iterating through 1e8 points to assign each one a probability takes forever (specifically because I need to make 1e8 calls to an interpolation function), and moreover is unnecessary. I can jump to a 100X100 grid without really losing precision. This is not the case for RVs and astrometry, where the likelihood of a given model depends on all of its orbital parameters.



def abs_Xmag(d_star, vmag, imaging_wavelength):
    """
    Calculate the absolute magnitude of the host star in the band nearest
    to the the wavelength of the observations.
    
    Arguments:
        d_star (float, au): Distance of system from Sun
        vmag (float, mag): Apparent V-band magnitude of host star
        imaging_wavelength (float, micrometers): Wavelength of imaging data
    
    returns:
        band_name (str): Name of band pass whose center is nearest to
                        the imaging wavelength
        host_abs_Xmag (float, mag): Absolute magnitude of host star in
                                    the band_name band
    """
    
    host_abs_vmag = vmag - 5*np.log10(d_star/pc_in_au/10) # absolute v-band magnitude. First convert au to pc
    
    # Find the central band that is closest to the imaging wavelength
    band_name = bands.iloc[np.argmin(abs(bands['wavelength']-imaging_wavelength))]['band']
    
    if band_name not in mamajek_table.columns:
        raise Exception('There is no data for the imaging band in the Mamajek table.\
                         Find imaging data in the V, R, I, J, H, or K bands')
    else:
        # Linearly interpolate between Mamajek entries to estimate Vmag ==> Xmag conversion
    
        f = interp1d(mamajek_table['V'], mamajek_table[band_name])
    
    host_abs_Xmag = f(host_abs_vmag)

    return band_name, host_abs_Xmag


def mag_to_mass(mags, masses, abs_Xmag_list):
    """
    Convert a list of absolute magnitudes in some band into companion masses.
    We do this by first interpolating between two corresponding lists of
    magnitudes and masses to create a function. Then apply this function to
    the list of absolute magnitudes to calculate corresponding masses.
    
    mags (list of floats, mag): List of reference absolute magnitudes in a 
                                given band
    masses (list of floats, M_jup): List of companion masses with the same 
                                   length as mags
    abs_Xmag_list (list of floats, mag): List of absolute magnitudes from the
                                        contrast curve.
    
    
    returns:
        mass_list (list of floats, M_jup): Interpolated masses corresponding to
                                          the magnitudes in abs_Xmag_list
    """
    
    # First, create an interpolated function to get mass from magnitude
    f = interp1d(mags, masses) # Function takes mag to mass
    
    mass_list = f(abs_Xmag_list)
    
    return mass_list
    


# if __name__ == "__main__":

# import matplotlib.pyplot as plt
# ####################################
# contrast_curve = pd.read_csv('data/TOI1174_832_sensitivity.dat',
#                               skiprows=29,
#                               delimiter=' ',
#                               header=None)
# new_header = ['ang_sep', 'delta_mag']
# contrast_curve.columns = new_header
# d_star = 94.9
# vmag = 10.96
# wavelength = 0.832
# ####################################
#
# my_interp = m_a_interp_fn(d_star, vmag, wavelength, contrast_curve)
#
# a_list = np.linspace(1, 90, 100)
#
# plt.plot(a_list, my_interp(a_list))
# plt.show()






















