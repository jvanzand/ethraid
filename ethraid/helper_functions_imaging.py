import os
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

from ethraid import _ROOT
import ethraid.compiled.helper_functions_general as hlp

# This module makes use of Table 5 of Pecaut & Mamajek (2013), available at https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
# It also uses Table 4 from Baraffe+03
# Both of these tables are included in the data/ directory
# Estimate that Ks (2.15μm) ~ K (2.2μm) and W1 (3.37μm) ~ L' (3.77μm) because Baraffe uses K and L'
mamajek_rename = {'K_s':'K', 'W1':'L_prime'}
mamajek_table = pd.read_csv(_ROOT+'/data/mamajek.csv').rename(columns=mamajek_rename)
baraffe_table = pd.read_csv(_ROOT+'/data/baraffe_table_4.csv')
bands= pd.read_csv(_ROOT+'/data/bands.csv')
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value


def imag_array(d_star, vmag, imag_wavelength, contrast_str, a_lim, m_lim, grid_num):
    """
    Convert a contrast curve in (angular_separation, Δmag) space
    into (separation, mass) space. This requires finding a
    correspondence between brightness and mass, which I do with
    either the Mamajek table or the Baraffe (03) paper, depending
    on the mass regime.
    
    Arguments:
        d_star (float, AU): Distance to the host star
        vmag (float, mag): Apparent V-band magnitude of host star
        imag_wavelength (float, μm): Wavelength of imaging data in contrast_curve
        contrast_str (str): Path to contrast curve with columns of 
                              'ang_sep' (arcseconds) and 'delta_mag' (mag)
        a_lim (tuple of floats, au): Semi-major axis limits to consider, 
                                     in the form (a_min, a_max).
        m_lim (tuple of floats, M_jup): Mass limits, (m_min, m_max).
        grid_num (int): Dimensions of (a,m) grid.
    
    returns:
        np_imag_array (2D array, dim=grid_numXgrid_num): Normalized probability
                                        array from imaging data. In this model,
                                        imaging probabilities are either 0 or
                                        "1" (actually 1 / the sum of the 
                                        array before normalization).
    """
    # If no imaging is provided, just return an array of 1s
    if vmag is None or imag_wavelength is None or contrast_str is None:
        print("helper_functions_imaging.imag_array: At least one imaging param is None. \n"
              "                                     post_imag will be a uniform array.")
        return np.ones((grid_num, grid_num))

    # Get band of observations and host star abs mag in that band
    band_name, host_abs_Xmag = abs_Xmag(d_star, vmag, imag_wavelength)

    contrast_curve = pd.read_csv(contrast_str)

    if not set(['ang_sep', 'delta_mag']).issubset(contrast_curve.columns):
        raise Exception("The dataframe must contain columns 'ang_sep' and 'delta_mag'")

    # Objective is sep (AU) vs. mass (M_jup)
    ############## 1: Get physical separation. Make sure to convert d_star from au to pc
    seps = contrast_curve['ang_sep']*(d_star/pc_in_au)
    ##############

    ############## 2: Get companion absolute mag at each Δmag
    # The absolute magnitudes of the companions are the delta_mag values plus the host star abs mag
    Xband_abs_mags = contrast_curve['delta_mag']+host_abs_Xmag

    ############## 3: Create interp_df, a csv relating companion absolute magnitude to companion mass
    # Start with Mamajek. We could cut super bright targets from the interpolation, but no need
    interp_df_mamajek = mamajek_table[['M_jup', band_name]]

    # It may be that Mamajek covers the whole contrast curve and we don't need Baraffe, but maybe not.
    # Find the dimmest mag in Mamajek and start using Baraffe mags after that.
    max_mag = interp_df_mamajek[band_name].max()

    # Just like we took even the brightest entries from Mamajek, take even the dimmest from Baraffe to be safe.
    interp_df_baraffe = baraffe_table.query("{}>{}".format(band_name, max_mag))[['M_jup', band_name]]

    # Concatenate the two dfs above
    interp_df = pd.concat([interp_df_mamajek, interp_df_baraffe]).sort_values(by=band_name)
    
    ############## 4: Use interpolation function to get companion masses from Xband_abs_mags
    companion_masses = mag_to_mass(interp_df[band_name], interp_df['M_jup'], Xband_abs_mags)
    ##############
    
    # Discrete correspondence between separation and mass. We get a continuous correspondence below
    a_m_contrast = pd.DataFrame({'M_jup':companion_masses, 'sep':seps}).reset_index(drop=True)

    
    # We might want to plot sma values greater than what's given in the contrast curve.
    # In that case, conservatively estimate that the curve becomes flat after the last sma value
    # (It would actually continue to drop to lower masses at larger separations, but increasingly slowly, so this is a fine approximation)
    # Similarly, what if we plot sma values BELOW the lowest contrast value?
    # Again, conservatively estimate the contrast becomes -inf, meaning that imaging rules out NO companions at those separations.
    last_a_ind = a_m_contrast[a_m_contrast['sep'] == a_m_contrast['sep'].max()].index # Find largest a
    last_m = a_m_contrast.iloc[last_a_ind]['M_jup'] # Find m corresponding to largest a
    
    # Mass values below bounds are +inf (bc -inf contrast), while those above bounds are last m value within bounds.
    a_m_interp_fn = interp1d(a_m_contrast['sep'], a_m_contrast['M_jup'], 
                             bounds_error=False, fill_value=(np.inf,last_m))
                             
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


# An important question I might ask later: why do I jump straight to “imag_array” (100X100) rather than first creating “imag_list” (len=1e8) and forming it into an array, as I do for RVs and astrometry? Answer: Basically, it takes too long and it's not needed. The constraints my model imposes on imaging are fairly rough: at a given separation, rule out all masses greater than X. Iterating through 1e8 points to assign each one a probability takes forever (specifically because I need to make 1e8 calls to an interpolation function), and moreover is unnecessary because I model imaging constraints as dependent ONLY on mass and sma. I can jump to a 100X100 grid without losing precision. This is not the case for RVs and astrometry, where the likelihood of a given model depends on all of its orbital parameters.

def abs_Xmag(d_star, vmag, imaging_wavelength):
    """
    Calculate the absolute magnitude of the host star in the band nearest
    to the the wavelength of the observations.
    
    Arguments:
        d_star (float, au): Distance of system from Sun
        vmag (float, mag): Apparent V-band magnitude of host star
        imaging_wavelength (float, micrometers): Wavelength of imaging data
    
    Returns:
        band_name (str): Name of band pass whose center is nearest to
                        the imaging wavelength
        host_abs_Xmag (float, mag): Absolute magnitude of host star in
                                    the band_name band
    """
    
    host_abs_vmag = vmag - 5*np.log10(d_star/pc_in_au/10) # absolute V-band magnitude. First convert au to pc
    
    # Find the central band that is closest to the imaging wavelength
    band_name = bands.iloc[np.argmin(abs(bands['wavelength']-imaging_wavelength))]['band']

    if band_name not in mamajek_table.columns:
        raise Exception("There is no data for the imaging band in the Mamajek table.\
                         Find imaging data in the V, R, I, J, H, K, or L' bands")
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
    
    Arguments:
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


def interp_fn(d_star, vmag, imag_wavelength, which='C2M'):
    """
    Helper function to interpolate between contrast and mass. Choose
    whether the output function takes contrasts and returns mass, or
    vice versa.
    
    Arguments:
        d_star (float, AU): Distance to the host star
        vmag (float, mag): Apparent V-band magnitude of host star
        imag_wavelength (float, μm): Wavelength of imaging data in contrast_curve
        which (str): Either 'C2M' or 'M2C' to choose whether contrast/mass
                     should be the input/output of the returned function.
    
    Returns:
        interp_fn (scipy function): Function taking contrast as an argument
                                    and returning companion mass, or vice versa.
    """
    
    # Get band of observations and host star abs mag in that band
    band_name, host_abs_Xmag = abs_Xmag(d_star, vmag, imag_wavelength)
    
    ############## Create interp_df, a csv relating companion absolute magnitude to companion mass
    # Start with Mamajek. We could cut super bright targets from the interpolation, but no need
    interp_df_mamajek = mamajek_table[['M_jup', band_name]]
    
    # It may be that Mamajek covers the whole contrast curve and we don't need Baraffe, but maybe not.
    # Find the dimmest mag in Mamajek and start using Baraffe mags after that. 
    max_mag = interp_df_mamajek[band_name].max()

    # Just like we took even the brightest entries from Mamajek, take even the dimmest from Baraffe to be safe.
    interp_df_baraffe = baraffe_table.query("{}>{}".format(band_name, max_mag))[['M_jup', band_name]]

    # Concatenate the two dfs above
    interp_df = pd.concat([interp_df_mamajek, interp_df_baraffe]).sort_values(by=band_name)[[band_name, 'M_jup']]
    interp_df['delta_mag'] = interp_df[band_name] - host_abs_Xmag # Contrast column
    
    if which=='C2M': # Function from contrast to mass
        interp_fn = interp1d(interp_df['delta_mag'], interp_df['M_jup'])
    
    elif which=='M2C': # Function from mass to contrast
        interp_fn = interp1d(interp_df['M_jup'], interp_df['delta_mag'])
        
    
    return interp_fn

















