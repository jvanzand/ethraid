import os
from ethraid import _ROOT, Ms2Mj, pc_in_au

# GENERAL PARAMS
# Number of orbital models to simulate
num_points = 1e6
# Dimension of grid over which model probabilities will be spread.
# Higher grid_num gives greater resolution, but fewer models per grid box.
# When using CLI, grid_num should be supplied at the command line rather than config file (default grid_num=100).
# Include grid_num in config file for API usage.
grid_num = 100
# Minimum and maximum semi-major axes to sample (AU)
min_a = 2
max_a = 1e2
# Minimum and maximum masses to sample (M_Jup)
min_m = 1
max_m = 1e3
# Eccentricity distribution for sampled orbits
e_dist = 'piecewise'


# STELLAR PARAMS
# Star name for file labeling. Need not be official.
star_name = '12572'
# Mass of star in Jupiter masses
m_star = 0.91*Ms2Mj
# Distance from Earth to star in AU
d_star = 65.9*pc_in_au


# RV PARAMS
# Whether to use RV data. Assign run_rv=False to omit RVs from the calculation entirely.
run_rv = True
# Linear RV trend term (m/s/day).
gammadot = -0.0599
# Error on gammadot
gammadot_err = 0.0037
# Quadratic RV curvature term (m/s/day/day)
gammaddot = 2.2e-6
# Error on gammaddot
gammaddot_err = 4.9e-6
# Epoch at which gammadot and gammaddot are measured. Typically about 1/2 way through the observing baseline.
rv_epoch = 2458991.236308


# ASTROMETRY PARAMS
# Whether to use astrometry data. Assign run_astro=False to omit astrometry from the calculation entirely.
run_astro = True
# Difference between the average Gaia proper motion and the position-based average proper motion between the Hipparcos and Gaia missions (milli-arcseconds/year)
# Set dmu/dmu_err to None to provide Hipparcos or Gaia ID instead
delta_mu = None
# Error on delta_mu
delta_mu_err = None
# Target Hipparcos identifier. Alternative to supplying delta_mu
hip_id = '9618'
# Target Gaia DR3 identifier. Alternative to supplying delta_mu
gaia_id = None


# IMAGING PARAMS
# Whether to use imaging data. Assign run_imag=False to omit imaging from the calculation entirely.
run_imag = True
# How to calculate imaging posterior. If 'exact', forward model companions as with RVs and astrometry.
# If 'approx', then for the imaging calculations only, approximate all orbits to be face-on and circular regardless of sampled parameters, and rule out any model with a mass/angular separation combo that was detectable by imaging.
imag_calc = 'exact'
# Host star visual magnitude. Used to estimate the magnitude at the imaging wavelength
vmag = 9.2
# Wavelength at which contrast curve was acquired (micrometers)
imag_wavelength = 2.2
# Path to contrast curve
contrast_str = os.path.join(_ROOT, 'data/clean_curves/TOI1471_Brgamma.csv')
# Epoch of imaging observations (BJD).
imag_epoch = 2458991.236308

# SAVE PARAMS
# Whether to save raw arrays (1d, unbinned), processed arrays (2d, binned), or both
save = ['raw', 'proc']
# Output directory. Destination of folder containing saved products
outdir = ''


# PLOTTING PARAMS
# Coordinates at which to plot a gold star. Usually corresponds to a known companion which could be the source of an observed trend.
scatter_plot = None