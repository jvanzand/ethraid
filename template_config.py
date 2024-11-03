## Template config file

from ethraid import Ms2Mj, pc_in_au
import ethraid.compiled.helper_functions_general as hlp

# GENERAL PARAMS
# Number of orbital models to simulate
num_points = 1e4

# Dimension of grid over which model probabilities will be spread.
# Higher grid_num gives greater resolution, but fewer models per grid box.
grid_num = 100
## Limits of parameter space to search. a = separation (AU), m = mass (M_Jup)
min_a = 2
min_m = 2

max_a = 1e2
max_m = 1e3

# Eccentricity distribution for sampled orbits
e_dist = 'piecewise'


# STELLAR PARAMS
# Star name for file labeling. Not used for catalog cross-matching
star_name = 'template_star'
# Mass of star in Jupiter masses. Ms2Mj converts M_Sun to M_Jup.
m_star = 0.5*Ms2Mj
# Distance from Earth to star in AU
d_star = 5*pc_in_au


# RV PARAMS
# Whether to use RV data. Assign run_rv=False to omit RVs from the calculation entirely.
run_rv = True
# Linear RV trend term (m/s/day).
gammadot = 0.1383
# Error on gammadot
gammadot_err = 0.0029
# Quadratic RV curvature term (m/s/day/day)
gammaddot = -7.44e-5
# Error on gammaddot
gammaddot_err = 2.7e-5
# Epoch at which gammadot and gammaddot are measured. Typically about 1/2 way through the observing baseline.
rv_epoch = 2458847.780463


# ASTROMETRY PARAMS
# Whether to use astrometry data. Assign run_astro=False to omit astrometry from the calculation entirely.
run_astro = True

# Difference between the average Gaia proper motion and the position-based average proper motion between the Hipparcos and Gaia missions (milli-arcseconds/year)
# Set dmu/dmu_err to None to provide Hipparcos or Gaia ID instead
delta_mu = 0.13
# Error on delta_mu
delta_mu_err = 0.03
# Target Hipparcos identifier. Alternative to supplying delta_mu
hip_id = None
# Target Gaia DR3 identifier. Alternative to supplying delta_mu
gaia_id = None


# IMAGING PARAMS
# Whether to use imaging data. Assign run_imag=False to omit imaging from the calculation entirely.
run_imag = True
# How to calculate imaging posterior. If 'exact', forward model companions as with RVs and astrometry.
# If 'approx', then for the imaging calculations only, approximate all orbits to be face-on and circular regardless of sampled parameters, and rule out any model with a mass/angular separation combo that was detectable by imaging.
imag_calc = 'exact'
# Host star visual magnitude. Used to estimate the magnitude at the imaging wavelength
vmag = 8.97
# Wavelength at which contrast curve was acquired (micrometers)
imag_wavelength = 2.2
# Choose (based on system age) which table from Baraffe+03 to use for brown dwarf cooling model
# Table 1-->0.1 Gyr, 2-->0.5 Gyr, 3-->1 Gyr, 4-->5 Gyr, 5-->10 Gyr
age_table = 4
# Path to contrast curve
contrast_str = os.path.join(_ROOT, 'data/test_K_band_altered.csv')
# Epoch of imaging observations (BJD).
imag_epoch = 2458795.5

# SAVE PARAMS
# Whether to save raw arrays (1d, unbinned), processed arrays (2d, binned), or both
save = ['raw', 'proc']
# Output directory. Destination of folder containing saved products
outdir = ''


# PLOTTING PARAMS
# Coordinates at which to plot a gold star (AU, M_Jup). Usually corresponds to a known companion which could be the source of an observed trend.
scatter_plot = [(5, 5)]