## Test config file 1. Runs RVs, astrometry (using delta_mu directly), and imaging (exact)
from ethraid import Ms2Mj, pc_in_au
import ethraid.compiled.helper_functions_general as hlp

# GENERAL PARAMS
# Number of orbital models to simulate
num_points = 1e4
# Dimension of grid over which model probabilities will be spread.
# Higher grid_num gives greater resolution, but fewer models per grid box.
# When using CLI, grid_num should be supplied at the command line rather than in the config file (default grid_num=100).
# Include grid_num in config file for API usage.
grid_num = 100
# Whether to save raw arrays (1d, unbinned), processed arrays (2d, binned), or both
save = ['raw', 'proc']
# Output directory. Destination of folder containing saved products
outdir = ''
## See below for min_a and min_m, calculated in terms of the RV data


# STELLAR PARAMS
# Star name for file labeling. Need not be official.
star_name = 'test1'
# Mass of star in Jupiter masses
m_star = 0.807*Ms2Mj
# Distance from Earth to star in AU
d_star = 53.8*pc_in_au


# RV PARAMS
# Whether to use RV data. Assign run_rv=False to omit RVs from the calculation entirely.
run_rv = True
# Linear RV trend term (m/s/day).
gammadot = 0.114
# Error on gammadot
gammadot_err = 0.006
# Quadratic RV curvature term (m/s/day/day)
gammaddot = -6e-5
# Error on gammaddot
gammaddot_err = 1.9e-5
# Epoch at which gammadot and gammaddot are measured. Typically about 1/2 way through the observing baseline.
rv_epoch = 2458847.780463


###################
## Calculate min a and min m based on RVs, so min_a, max_a, min_m, and max_m must be defined after the RV params
rv_baseline=500
min_per=1000
min_a_and_m = hlp.min_a_and_m(gammadot, gammaddot, rv_baseline, min_per, m_star)

min_a = min_a_and_m[0]
min_m = min_a_and_m[1]

max_a = 1e2
max_m = 1e3
###################


# ASTROMETRY PARAMS
# Whether to use astrometry data. Assign run_astro=False to omit astrometry from the calculation entirely.
run_astro = True
# Difference between the average Gaia proper motion and the position-based average proper motion between the Hipparcos and Gaia missions (milli-arcseconds/year)
# Set dmu/dmu_err to None to provide Hipparcos or Gaia ID instead
delta_mu = 0.12767382507786398
# Error on delta_mu
delta_mu_err = 0.034199052901953214
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
# Path to contrast curve
contrast_str = 'ethraid/data/test_K_band.csv'
# Epoch of imaging observations (BJD).
imag_epoch = 24593300


# PLOTTING PARAMS
# Coordinates at which to plot a gold star. Usually corresponds to a known companion which could be the source of an observed trend.
scatter_plot = [3.8, 3.0]