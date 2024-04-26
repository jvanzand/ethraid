from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

cimport numpy as np
from libc.math cimport sin, cos, sqrt
from ethraid.compiled._kepler import kepler_single

from ethraid import _ROOT, hip_times
import ethraid.compiled.helper_functions_general as hlp

cdef double two_pi, pc_in_au

two_pi = 6.283185307179586
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value

# This module makes use of Table 5 of Pecaut & Mamajek (2013), available at https://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt
# It also uses Table 4 from Baraffe+03
# Both of these tables are included in the data/ directory
# Estimate that Ks (2.15μm) ~ K (2.2μm) and W1 (3.37μm) ~ L' (3.77μm) because Baraffe uses K and L'
mamajek_rename = {'K_s':'K', 'W1':'L_prime'}
mamajek_table = pd.read_csv(_ROOT+'/data/mamajek.csv')\
                  .rename(columns=mamajek_rename)
#baraffe_table = pd.read_csv(_ROOT+'/data/baraffe_table4.csv')
bands= pd.read_csv(_ROOT+'/data/bands.csv')



def imag_list(double [:] a_list, double [:] m_list, double [:] e_list, 
              double [:] i_list, double [:] om_list, 
              double [:] M_anom_0_list, double [:] per_list, 
              double m_star, double d_star, double vmag,
              double imag_wavelength, int age_table, double imag_epoch, str contrast_str):
    """
    Calculates the likelihood of the imaging data conditioned 
    on each of a list of orbital models, resulting in a list of 
    likelihoods. All input lists must have the same length, which 
    is also the length of the output, log_lik_list.
    
    NOTE that unlike for RVs and astrometry, each log-likelihood is either
    0 or -np.inf (likelihood = 1 or 0, respectively), approximating that 
    with respect to the imaging data, a given model has a 100% or 0% chance 
    of being detected.
    
    Arguments:
        a_list (list of floats, AU): Semi-major axes
        m_list (list of floats, M_jup): Companion masses
        e_list (list of floats): Orbital eccentricities
        i_list (list of floats, radians): Orbital inclinations
        om_list (list of floats, radians): Arguments of periastron
        M_anom_0_list (list of floats, radians): Mean anomalies at the 
                                                 beginning of the Hipparcos 
                                                 mission (1989.85)
        per_list (list of floats, days): Orbital periods
        m_star (float, M_jup): Host star mass
        d_star (float, AU): Distance of system from Earth
        vmag (float): Apparent V-band magnitude of host star
        imag_wavelength (float, micrometers): Wavelength of imaging data
        age_table (int): Integer 1-5, indicating which BD cooling model to use
                         based on age of system.
                         1-->0.1 Gyr, 2-->0.5 Gyr, 3-->1 Gyr, 4-->5 Gyr, 5-->10 Gyr
        imag_epoch (float, BJD): Epoch at which imaging was acquired
        contrast_str (str): Path to contrast curve file
    
    Returns:
        log_lik_list (list of floats): List of likelihoods corresponding
                                       to the input models. 
    """
    
    cdef int num_points, j
    cdef double a, m, e, i, om, M_anom_0, per
    
    num_points = a_list.shape[0]
    cdef np.ndarray[double, ndim=1] ang_sep_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    max_dmag_list = np.ndarray(shape=(num_points,), dtype=np.float64),\
                                    model_dmag_list = np.ndarray(shape=(num_points,), dtype=np.float64)                               
    
    
    print('Running imaging models')
    for j in tqdm(range(num_points)):
    
        a = a_list[j]
        m = m_list[j]
        e = e_list[j]
        i = i_list[j]
        om = om_list[j]
        M_anom_0 = M_anom_0_list[j]
        per = per_list[j]
        
        # Angular separation in arcseconds added to list
        ang_sep_list[j] = ang_sep(a, m, e, i, om, M_anom_0, per, 
                                  m_star, d_star, imag_epoch)
                                  
    ##########
    # This is the "data"
    # Interpolate a contrast curve to get a function that maps angular separation to delta_mag contrasts
    angsep_to_dmag = interp_fn(d_star, vmag, imag_wavelength, age_table, 
                               contrast_str=contrast_str, which='A2C')
    
    # This is the "model"
    # Next function is to map model companion masses to delta_mag contrasts based on stellar/Brown Dwarf mass-luminosity models from Pecaut/Mamajek2013 and Baraffe+03.
    # Must set fill_value for mass_to_dmag because sampled masses go well below masses listed in Brown Dwarf cooling model tables (m_min~2 M_J). For objects smaller than that, make contrast very high (ie, we cannot image 2 M_J objects).
    mass_to_dmag = interp_fn(d_star, vmag, imag_wavelength, age_table, which='M2C', fill_value=(10000, 0))
    #########
    
    # The "data" sets the dimmest detectable object at the model seps
    max_dmag_list = angsep_to_dmag(ang_sep_list)
    # The "model" gives the brightness of the modeled companion
    model_dmag_list = mass_to_dmag(m_list)
    
    log_lik_list = model_dmag_list > max_dmag_list # lik_list is a list of booleans. If target is too dim, False
    log_lik_list = log_lik_list.astype(float) # Convert booleans to 1s and 0s
    log_lik_list[log_lik_list==0] = -np.inf # Convert lik to log-lik
    log_lik_list[log_lik_list==1] = 0 # Convert lik to log-lik
    
    return log_lik_list


def ang_sep(double a, double m, double e, double i, double om, double M_anom_0, 
            double per, double m_star, double d_star, double imag_epoch):
    """
    Computes the projected angular separation of a companion on a 
    model orbit. This function is analogous to dmu() in the astrometry
    module.
    
    Arguments:
        a (float, AU): Semi-major axis
        m (float, M_jup): Companion mass
        e (float): Orbital eccentricity
        i (float, radians): Orbital inclination
        om (float, radians): Argument of periastron
        M_anom_0 (float, radians): Mean anomaly at the beginning
                                   of the Hipparcos mission (1989.85)
        per (float, days): Orbital period
        m_star (float, M_jup): Host star mass
        d_star (float, AU): Distance of system from Earth
        imag_epoch (float, BJD): Date at which imaging data was
                                 acquired
    
    Returns:
        ang_sep (float, arcsec): Modeled angular separation
                                 at the imaging epoch.
        
    """
    cdef double mass_ratio, au_2_mas, mean_motion, M, E, ang_sep
    
    cdef double [:,:] rot_mtrx # This makes rot_mtrx a memview
    rot_mtrx = np.zeros((3,3),dtype=np.float64)
    
    # vec holds various values throughout ang_sep(). After each value
    # has served its purpose, it is overwritten so that only one
    # vector needs to be allocated, saving time.
    cdef double vec_list[3]
    cdef double [:] vec = vec_list
    
    
    mass_ratio = m/(m_star + m)
    d_star = d_star/206264.80624548031 # Divide by (c.pc.cgs/c.au.cgs).value to get units of pc
    au_2_as = 1/d_star # Conversion factor btwn au and arcseconds, with d_star converted to pc

    mean_motion = two_pi/per
    rot_matrix(i, om, rot_mtrx)
    
    # Advance initial mean anomaly from the beginning of the Hipparcos mission to the imaging epoch
    M = mean_motion*(imag_epoch-hip_times[0]) + M_anom_0
    
    E = kepler_single(M%two_pi, e)
    
    # Equations 2.41 in Murray & Dermott
    # vec points from the star to the companion (note the pre-factor, a, is the full semi-major axis)
    vec[0] = a*(cos(E)-e)
    vec[1] = a*sqrt(1-e*e)*sin(E)
    vec[2] = 0
    
    # mat_mul() returns a vector describing stellar position in the observer's frame (in the sky plane, with the z-axis pointing to the observer)
    # mat_mul() takes vec as input, overwrites it, and returns it as output. This saves time.
    mat_mul(rot_mtrx, vec, vec)
    
    # Calculate magnitude of only the x and y components of vec because we want the projected separation on the sky.
    ang_sep = sqrt((vec[0]*au_2_as)**2 + (vec[1]*au_2_as)**2)
    
    return ang_sep
    
    
    
cdef void rot_matrix(double i, double om, double [:,:] rot_mtrx):
    """
    This is P2*P1 from Murray & Dermott Eqs. 2.119. We omit P3 (2.120) because
    the longitude of the ascending node (Omega) can be set arbitrarily to 0 for
    our purposes, saving some time.

    This function doesn't return anything. Instead, declare a matrix in your 
    function and this will update it, saving lots of time by not allocating memory 
    to and returning a matrix.

    Arguments:
        i (float, radians): orbital inclination (0 = face-on)
        om (float, radians): argument of periastron
        Om (float, radians): longitude of the ascending node
        rot_mtrx (3x3 array of zeros): Initial array, modified in place

    Returns:
        None (but populates rot_mtrx with values)
    """
    cdef double sin_om, sin_i, cos_om, cos_i

    sin_om = sin(om)
    sin_i  = sin(i)
    cos_om = cos(om)
    cos_i  = cos(i)

    rot_mtrx[0][0] = cos_om
    rot_mtrx[0][1] = -sin_om
    rot_mtrx[0][2] = 0

    rot_mtrx[1][0] = cos_i*sin_om
    rot_mtrx[1][1] = cos_i*cos_om
    rot_mtrx[1][2] = -sin_i

    rot_mtrx[2][0] = sin_i*sin_om
    rot_mtrx[2][1] = sin_i*cos_om
    rot_mtrx[2][2] = cos_i

    # Do not return rot_mtrx
        

# THIS version of the function will let me input the same vector as in_vec and out_vec.
# I want to do this because it will allow me to allocate fewer arrays in the ang_sep 
# function, which offers speedups.
cpdef void mat_mul(double [:,:] mat, double [:] in_vec, double [:] out_vec):
    """
    This function is written specifically to matrix multiply rot_mtrx (3x3)
    with vec (3x1) in the dmu() function. Saves time by receiving its 
    "output" (out_vec) as an argument and modifying it in place without
    returning anything. I made this function even more ad-hoc by removing
    the for-loop over vector elements. This lets the user pass the same
    vector as both the input and the output, which saves time by only
    allocating one array.

    Arguments:
        mat_mul (3x3 array of floats): Rotation matrix from the rot_matrix
                                       function
        in_vec (3x1 array of floats): Position vector of star [x, y, z]
    
    Returns:
        None (but populates out_vec)
    """

    cdef double m00, m01, m02,\
                m10, m11, m12,\
                m20, m21, m22,\
                iv0, iv1, iv2,\
                ov0, ov1, ov2

    m00 = mat[0][0]
    m01 = mat[0][1]
    m02 = mat[0][2]
    m10 = mat[1][0]
    m11 = mat[1][1]
    m12 = mat[1][2]
    m20 = mat[2][0]
    m21 = mat[2][1]
    m22 = mat[2][2]

    iv0 = in_vec[0]
    iv1 = in_vec[1]
    iv2 = in_vec[2]

    ov0 = m00*iv0 + m01*iv1 + m02*iv2
    ov1 = m10*iv0 + m11*iv1 + m12*iv2
    ov2 = m20*iv0 + m21*iv1 + m22*iv2

    out_vec[0] = ov0
    out_vec[1] = ov1
    out_vec[2] = ov2
        
    # Do not return out_vec
    
    
def interp_fn(d_star, vmag, imag_wavelength, age_table, contrast_str=None, which='C2M', fill_value=None):
    """
    Helper function to interpolate between Δmag contrast and mass. 
    Choose whether the output function takes contrasts and returns 
    mass, or vice versa.

    Arguments:
        d_star (float, AU): Distance to the host star
        vmag (float, mag): Apparent V-band magnitude of host star
        imag_wavelength (float, μm): Wavelength of imaging data in contrast_curve
        age_table (int): Integer from 1 to 5, indicating which BD cooling model to use
                         based on age of system. Correspond to tables 1-5 in Baraffe+03.
                         1-->0.1 Gyr, 2-->0.5 Gyr, 3-->1 Gyr, 4-->5 Gyr, 5-->10 Gyr
        contrast_str (str): Path to a csv file with columns "ang_sep" and "delta_mag"
        which (str): One of 'C2M', 'M2C', 'A2M', or 'A2C' to choose a function with
                     arguments/outputs of contrast/mass, mass/contrast,
                     angular separation/mass, or angular separation/contrast
        fill_value (float or tuple): Value to use when the interpolation function
                                     is called with an input outside of the original
                                     interpolation range. If float, that value will
                                     be used for any argument outside the range. If
                                     2-tuple, first value is used for arguments below
                                     interpolation range, and second for args above.
                                     For example, if interp_fn is to accept angular
                                     separation and output companion mass, fill_value
                                     could be np.inf, to reflect that above/below the
                                     known contrast curve, you cannot rule out any
                                     companion masses.

    Returns:
        interp_fn (scipy function): Function taking contrast as an argument
                                    and returning companion mass, or vice versa.
                                    If contrast_str is also provided, a third
                                    option 'A2M' takes angular separation and
                                    returns companion mass.
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
    baraffe_table = pd.read_csv(_ROOT+'/data/baraffe_table{}.csv'.format(age_table)) # Load correct table
    interp_df_baraffe = baraffe_table.query("{}>{}".format(band_name, max_mag))[['M_jup', band_name]]

    # Concatenate the two dfs above
    interp_df = pd.concat([interp_df_mamajek, interp_df_baraffe]).sort_values(by=band_name)[[band_name, 'M_jup']]
    interp_df['delta_mag'] = interp_df[band_name] - host_abs_Xmag # Contrast column

    if which=='C2M': # Function from contrast to mass
        interp_fn = interp1d(interp_df['delta_mag'], interp_df['M_jup'], 
                             bounds_error=False, fill_value=fill_value)

    elif which=='M2C': # Function from mass to contrast
        interp_fn = interp1d(interp_df['M_jup'], interp_df['delta_mag'],
                             bounds_error=False, fill_value=fill_value)
    
    elif which=='A2M': # Function from angular separation to mass
        if contrast_str is None:
            raise Exception('helper_functions_imaging.interp_fn: \n'
                            '                                   contrast_str required to generate A2M function')
        contrast_curve = pd.read_csv(contrast_str) # First get your angsep/dmag curve
        
        try:
            dmag_to_mass = interp1d(interp_df['delta_mag'], interp_df['M_jup'], bounds_error=False) # Create dmag->mass fn using Mamajek/Baraffe
        except ValueError:
            print("helper_functions_imaging.interp_fn: \n"
                  "                                   Value error while creating dmag->mass function.\n"
                  "                                   Try truncating your contrast curve to omit very \n"
                  "                                   small or large Δmag values.")
        
        contrast_masses = dmag_to_mass(contrast_curve['delta_mag']) # Use dmag->mass fn to find masses for each angsep in your contrast curve
        
        interp_fn = interp1d(contrast_curve['ang_sep'], contrast_masses,
                             bounds_error=False, fill_value=fill_value) # Finally, use angseps and corresponding masses to define an angsep->mass fn
    
    elif which=='A2C': # Very basic function that just interpolates contrast curve
        if contrast_str is None:
            raise Exception('helper_functions_imaging.interp_fn: \n'
                            '                                   contrast_str required to generate A2M function')
        contrast_curve = pd.read_csv(contrast_str) # First get your angsep/dmag curve
        fill_value = (0, contrast_curve['delta_mag'][-1]) # Close sep: 0 contrast; distant sep: equal to last contrast
        
        interp_fn = interp1d(contrast_curve['ang_sep'], contrast_curve['delta_mag'], 
                            bounds_error=False, fill_value=fill_value)
        
    return interp_fn

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
                         Please provide imaging data in the V, R, I, J, H, K, or L' bands")
    else:
        try:
            # Linearly interpolate between Mamajek entries to estimate Vmag ==> Xmag conversion
            f = interp1d(mamajek_table['V'], mamajek_table[band_name])
        except ValueError:
            print('helper_functions_imaging.abs_Xmag: Host Vmag must be between -3 and 19.4')

    host_abs_Xmag = f(host_abs_vmag)

    return band_name, host_abs_Xmag


def imag_array(d_star, vmag, imag_wavelength, age_table, contrast_str, a_lim, m_lim, grid_num):
    """
    Convert a contrast curve in (angular_separation, Δmag) space
    into (separation, mass) space. This requires finding a
    correspondence between brightness and mass, which is done with
    either the Mamajek table or the Baraffe (03) paper, depending
    on the mass regime.
    
    NOTE: This function makes the approximation that angular separation
          depends solely on semi-major axis (ie, assumes face-on orbits
          with e=0). For a full treatment, use the imag_list function.

    Arguments:
        d_star (float, AU): Distance to the host star
        vmag (float, mag): Apparent V-band magnitude of host star
        imag_wavelength (float, μm): Wavelength of imaging data in contrast_curve
        age_table (int): Integer 1-5, indicating which BD cooling model to use
                         based on age of system.
                         1-->0.1 Gyr, 2-->0.5 Gyr, 3-->1 Gyr, 4-->5 Gyr, 5-->10 Gyr
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

    contrast_curve = pd.read_csv(contrast_str)


    if not set(['ang_sep', 'delta_mag']).issubset(contrast_curve.columns):
        raise Exception("The dataframe must contain columns 'ang_sep' and 'delta_mag'")

    seps = contrast_curve['ang_sep']*(d_star/pc_in_au)*two_pi/8 # π/4 correction factor to avg over i and E
    dmags = contrast_curve['delta_mag']

    ## Create an interpolation fn that takes delta mag contrast to companion mass
    ## Fill value: if dmag is extremely small, mass is large. If too large, then mass -> 0
    dmag_to_mass = interp_fn(d_star, vmag, imag_wavelength, age_table, which='C2M', fill_value=(np.inf,0))
    companion_masses = dmag_to_mass(dmags)

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


    # Note that this is a coarse approximation. The "mass" over the whole ith grid square is just the conversion of the i+0.5th index to mass (and same for sma). Use i+0.5 to get to the middle of the grid square. i would be the far left (or bottom).
    for i in range(grid_num): # on mass axis
        m = hlp.index2value(i+0.5, (0, grid_num), m_lim)
        for j in range(grid_num): # on a axis
            a = hlp.index2value(j+0.5, (0, grid_num), a_lim)
        
            max_m = a_m_interp_fn(a) # max_m is the most massive object that could still be found
        
            # If the model mass is more than that, it has ~0 chance of existing, bc imaging didn't find it
            if m > max_m:
                imag_array[i, j] = 0
            
    np_imag_array = np.array(imag_array) / np.array(imag_array).sum()

    return np_imag_array
    
    
    