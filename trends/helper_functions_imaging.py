import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


baraffe_table = pd.read_csv('data/baraffe_table_4.csv')
mamajek_table = pd.read_csv('data/mamajek.csv')
bands = pd.read_csv('data/bands.csv')


def m_a_contrast(dstar, vmag, wavelength, contrast_curve):
    """
    Convert a contrast curve in (angular_separation, Δmag) space
    into (separation, mass) space. This requires finding a
    correspondence between brightness and mass, which I do with
    either the Mamajek table or the Baraffe (03) paper, depending
    on the mass regime.
    
    dstar (float, pc): Distance to the host star
    vmag (float, mag): Apparent V-band magnitude of host star
    wavelength (float, μm): Wavelength of imaging data in contrast_curve
    contrast_curve (dataframe, columns of 'ang_sep' (arcseconds) and 
                    'delta_mag' (mag)): Ordered pairs of angular
                                        separation and Δmag.
    
    returns:
        m_a_contrast (dataframe, columns of 'Msun' and 'sep' (AU))
    """
    
    band_name, host_abs_Xmag = abs_Xmag(dstar, vmag, wavelength)
    
    seps = contrast_curve['ang_sep']*dstar
    Xband_abs_mags = contrast_curve['delta_mag']+host_abs_Xmag
    
    # Start with Mamajek. We could cut super bright targets from the interpolation, but no need
    interp_df_mamajek = mamajek_table[['Msun', band_name]]
    
    # It may be that Mamajek covers the whole contrast curve and we don't need Baraffe, but maybe not.
    # Find the dimmest mag in Mamajek and start selecting Baraffe mags after that. 
    max_mag = interp_df_mamajek[band_name].max()

    # Just like we took super bright entries from Mamajek, take even the dimmest from Baraffe to be safe.
    interp_df_baraffe = baraffe_table.query("{}>{}".format(band_name, max_mag))[['Msun', band_name]]

    # Concatenate the two dfs above
    interp_df = pd.concat([interp_df_mamajek, interp_df_baraffe]).sort_values(by=band_name)
    
    
    companion_masses = mag_to_mass(interp_df[band_name], interp_df['Msun'], Xband_abs_mags)
    
    m_a_contrast = pd.DataFrame({'Msun':companion_masses, 'sep':seps})
    
    return m_a_contrast


def abs_Xmag(dstar, vmag, imaging_wavelength):
    """
    Calculate the absolute magnitude of the host star in the band nearest
    to the the wavelength of the observations.
    
    vmag (float, mag): Apparent V-band magnitude of host star
    imaging_wavelength (float, micrometers): Wavelength of imaging data
    
    returns:
        band_name (str): Name of band pass whose center is nearest to
                        the imaging wavelength
        host_abs_Xmag (float, mag): Absolute magnitude of host star in
                                    the band_name band
    """
    
    host_abs_vmag = vmag - 5*np.log10(dstar/10) # absolute v-band magnitude
    
    # Find the central band that is closest to the imaging wavelength
    band_name = bands.iloc[np.argmin(abs(bands['wavelength']-imaging_wavelength))]['band']
    
    if band_name not in mamajek_table.columns:
        raise Exception('There is no data for the imaging band in the Mamajek table.\
                         Find imaging data in the V, R, I, J, H, or K_s bands')
    else:
        # Linearly interpolate between Mamajek entries to estimate Vmag ==> Xmag conversion
        interp_list = mamajek_table[band_name]
    
    f = interp1d(mamajek_table['V'], interp_list)
    
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
    masses (list of floats, Msun): List of companion masses with the same 
                                   length as mags
    abs_Xmag_list (list of floats, mag): List of absolute magnitudes from the
                                        contrast curve.
    
    
    returns:
        mass_list (list of floats, Msun): Interpolated masses corresponding to
                                          the magnitudes in abs_Xmag_list
    """
    
    # First, create an interpolated function to get mass from magnitude
    f = interp1d(mags, masses) # Function takes mag to mass
    
    mass_list = f(abs_Xmag_list)
    
    return mass_list



if __name__ == "__main__":
    
    ####################################
    contrast_curve = pd.read_csv('data/TOI1174_832_sensitivity.dat', 
                                  skiprows=29, 
                                  delimiter=' ', 
                                  header=None)
    new_header = ['ang_sep', 'delta_mag']
    contrast_curve.columns = new_header
    dstar = 94.9
    vmag = 10.96
    wavelength = 0.832
    ####################################

    m_a_contrast_df = m_a_contrast(dstar, vmag, wavelength, contrast_curve)

    plt.plot(m_a_contrast_df['sep'], m_a_contrast_df['Msun'])
    plt.show()






















