#[constants]
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
Ms2Mj = M_sun/M_jup

# general params
num_points = 1e6
grid_num = 100
save = ['proc']
outdir = ''
verbose = False
min_m = 1
min_a = 1


#stellar_params
star_name = 'default'
m_star = 1.0*Ms2Mj
d_star = 10.0*pc_in_au


#rv_params
run_rv = False
gammadot = 0
gammadot_err = 1e8
gammaddot = 0
gammaddot_err = 1e8
rv_baseline = 500.
rv_epoch = 2459000.


#astrometry_params
run_astro = False
# Set dmu/dmu_err to None to provide Hipparcos or Gaia ID instead
delta_mu = 0
delta_mu_err = 1e8
hip_id = None
gaia_id = None


#imaging_params
run_imag = False
imag_calc = 'exact'
vmag = 10
imag_wavelength = 2.2
contrast_str = None
imag_epoch = 2459000


#plotting_params
scatter_plot = None