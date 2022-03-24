import pandas as pd

M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
Ms2Mj = M_sun/M_jup

# params_star = (star_name, m_star, distance(AU), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), rv_epoch, delta_mu, delta_mu_err, 
#               vmag=None, imag_wavelength=None, contrast_curve=None, scatter_tuple=None)


############## Syetems for papers #########################
params_191939 = ('191939', 0.807*Ms2Mj, 58.3*pc_in_au, 0.1116, 0.0037, -3.44e-5, 5.1e-6, 
                778.855, 2459192.641, 0.12767382507786398, 0.034199052901953214)

# Potential interesting DG paper.
params_12572 = ('12572', 0.91*Ms2Mj, 65.9*pc_in_au, -0.0549, 0.0049, -1.18e-5, 8.9e-6,
                980, 2459281, 0.0748781, 0.0451)#, 198, 0.085*Ms2Mj)
###########################################################

params_191939_old = ('191939_old', 0.807*Ms2Mj, 58.3*pc_in_au, 0.114, 0.006, -6e-5, 1.9e-5, 
                430.2527364352718, 2458847.780463, 0.12767382507786398, 0.034199052901953214,
                8.97, 0.832, pd.read_csv('data/191939_832_clean.csv'))
             
# Joey's target. Baseline and epoch taken from time series on Slack. It's a Hip target, but not in HGCA.
# Astrometry calculated manually from astro values in Table 3 of https://arxiv.org/pdf/1906.02058.pdf
params_hd206893 = ('HD206893', 1.32*Ms2Mj, 40.77*pc_in_au, -70/365.25, 31/365.25, 0, 1e8,
                    565.43, 2457949.229, 0.50, 0.10)

# CLS star with known companion from Lea's table 5.
params_hd182488= ('HD182488', 0.96*Ms2Mj, 3.196e6, -0.005390, 0.000345, -9.5e-7, 9e-8, 
                7951.02, 2454980, 1.0397, 0.0261, 20.97, 0.04*Ms2Mj,)
                
# CLS star with known companion from Lea's table 5.
params_hd201091= ('HD201091', 0.64*Ms2Mj, 7.182422e5, -0.007369, 0.000355, -3.22e-7, 5.4e-8, 
                11886.71, 2452901, 4.763910, 0.240117, 82.71, 0.59*Ms2Mj)
                
# DG star
params_t001438= ('T001438', 0.86*Ms2Mj, 111*pc_in_au, 0.079, 0.017, 3e-5, 1.2e-4, 
                433, 2459330, 0, 'a')

# CLS star with known companion.
params_hd186408= ('HD186408', 1065.317938, 4459779.594497, -0.004586, 0.000198, 0, 1e8, 
                12009.069844, 2452962.453122, 0.190556, 0.040519, 840.84, 0.99*Ms2Mj)

# CLS star with known companion from Lea's table 5.
params_hd131156= ('HD131156', 935.993904, 1381916.161366, 0.065991, 0.001904, 0, 1e8, 
                11998.216066, 2452956.829933, 19.315683, 0.175772, 35, 0.65*Ms2Mj)

# CLS star with known companion from Lea's table 5.
params_hd40397= ('HD40397', 925.627738, 4785726.363004, -0.029115, 0.000287, 0, 1e8, 
                8042.117923, 2454859.828823, 2.058092, 0.03366, 60.84, 0.27*Ms2Mj)

params_t001174= ('T001174', 0.82*Ms2Mj, 94.9*pc_in_au, -0.116, 0.019, 4.7e-5, 7.2e-5, 
                748, 2459209, -1, 0,
                10.96, 0.832, pd.read_csv('data/TOI1174_832_clean.csv'))

# # Test for 0 trend
# params_TEST = ('TEST', 0.807*Ms2Mj, 58.3*pc_in_au, 0.114, 0.006, -6e-5, 1.9e-5,
#                 430.2527364352718, 2458847.780463, 0.12767382507786398, 0.034199052901953214)

################ Systems for Keck proposal ##############
# No trend
params_t001410= ('T001410', 0.77*Ms2Mj, 72.9*pc_in_au, -0.002, 0.0034, 0, 1e10, 
                780, 2459162, 0, 'no_astro_sdf')
                
params_hd156141= ('HD156141', 1.03*Ms2Mj, 73*pc_in_au, 0.0596, 0.0075, -1e-5, 1.5e-5, 
                651, 2459345, 0, 'no_astro_sdf')


################ Validation Systems #################
         
# GL758, an example star in Tim Brandt's Orvara code. Using this to compare results.
params_gl758 = ('GL758', 0.95*Ms2Mj, 15.5*pc_in_au, -0.00633, 0.00025, 0, 1e8, #-8.19e-7, 0.67e-7,
                8413.010, 2454995.123, 1.0397, 0.0261)
                
# HIP67246, a confirmed giant published by Sarah Blunt. I truncated the timeseries to get gdot/gddot as a test.
# gdot/gddot fit with scipy.optimize.curve_fit bc rvsearch gave weird results
params_hip67246 = (1.10*Ms2Mj, 30.6*pc_in_au, 0.01, 0.002, 4.53e-5, 0.72e-5,
                1126.94, 2457685.312, 0.135, 0.035)
         
# HIP63510 aka Ross 458, a binary star system (primary is 0.6 M_sun, secondary is 80 M_J.)
# https://arxiv.org/pdf/1103.3544.pdf
params_hip63510 = (0.6*Ms2Mj, 11.4*pc_in_au, 1.53, 3*0.153, -1.9e-4, 1.9e-5,
                2336.7, 2455715.384, 14.465, 0.23)
                

# Synthetic: m = 32, a = 10
params_synth = ('synthetic', 1249.7, 15267565.22, -0.0861595, 0.008, -8.58e-6, 10000000,
                750.3, 2454737.94, 0.08233, 0.0465)

# # # Synthetic: a = 17, m = 47, i = pi/2, tiny astro signal
# params_synth = (1*Ms2Mj, 10*pc_in_au, -0.02944, 0.02944/5, 1.88e-5, 1.88e-5/5,
#                 300, 2459000, 0.59, 0.59/5)
                
# # # Synthetic: a = 17, m = 47, i = pi/3, larger astro signal
# params_synth = (1*Ms2Mj, 10*pc_in_au, -0.0255, 0.0231/5, 1.63e-5, 1.63e-5/5,
#                 300, 2459000, 3.282, 3.283/5)
                
                
                
                
                
                
                