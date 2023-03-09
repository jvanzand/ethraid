M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
Ms2Mj = M_sun/M_jup

# params_star = (star_name, m_star, distance(AU), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), rv_epoch, delta_mu, delta_mu_err, 
#               vmag=None, imag_wavelength=None, contrast_str=None, scatter_tuple=[sma, mass])

############## Systems for papers #########################

params_191939 = ('191939',  0.807*Ms2Mj, 53.8*pc_in_au, 0.1261, 0.0031, -5.92e-5, 3.4e-6, 
                980, 2459330, 0.12767382507786398, 0.034199052901953214,
                8.97, 2.2, 'ethraid/data/clean_curves/TOI1339_K_band.csv', [3.825, 3.0])

# # From 191939 paper (~mid-2021)
params_191939_old = ('191939_old', 0.807*Ms2Mj, 53.8*pc_in_au, 0.114, 0.006, -6e-5, 1.9e-5,
                430.2527364352718, 2458847.780463, 0.12767382507786398, 0.034199052901953214,
                8.97, 2.2, 'ethraid/data/clean_curves/TOI1339_K_band.csv')

## 6 DG trend systems for 2023B proposal
params_t001174 = ('T001174', 858.61968, 19575475.289, -0.077994, 0.008711, 1.4103e-4, 2.8656280365e-5, 
                1061.139983, 2459449.7717, None, None,
                10.96, 0.832, 'ethraid/data/clean_curves/TOI1174_832.csv')
params_t001246 = ('T001246', 913.8698, 35118549.3008, -0.0144221, 0.0047521, 2.296341e-5, 5.39474e-6, 
                923.661, 2459412.968705, None, None,
                11.63, 2.2, 'ethraid/data/clean_curves/TOI1246_K_band.csv')
params_t001438 = ('T001438', 905.3817, 22869884.545, 0.08318, 0.0106, -8.19e-6, 2.4323e-5, 
                755.914, 2459436.956, None, None,
                10.96, 0.832, 'ethraid/data/clean_curves/TOI1438_832.csv')
params_t001710 = ('T001710', 1065.65, 16795827.9265, -0.0138144566, 0.004381451657, 9.3825382e-8, 1.7494e-5,
                890.330352, 2459590.713166, None, None, 
                9.66, 0.832, 'ethraid/data/clean_curves/TOI1710_832.csv')
params_156141 = ('156141', 1074.86435, 15067296.662, 0.07555, 0.00441, -4.364061e-5, 5.4814e-6,
                897.8, 2459471.3113, None, None, 
                8.86, 2.2, 'ethraid/data/clean_curves/TOI1742_Brgamma.csv')
params_93963 = ('93963', 1132.1293873, 17048459.38787, -0.0389811, 0.007586, -1.43945635e-4, 2.7222588113e-5,
                922.304977000225, 2459607.259834, None, None,
                9.034, 0.832, 'ethraid/data/clean_curves/TOI1797_832.csv')


# Some DG trend systems
params_t001279 = ('T001279', 0.85*Ms2Mj, 107.4*pc_in_au, -0.0094, 0.003, 0, 1e8,
                817, 2459210, None, None, 
                10.71, 0.832, 'ethraid/data/clean_curves/T001279_832_speckle.csv')
params_t001422 = ('T001422', 1031.990824, 32143494.817747, 0.0128, 0.0035, 0, 1e8,
                828.887762, 2459217.272573, None, None, 
                10.62, 2.2, 'ethraid/data/clean_curves/T001422_2200_AO.csv')
params_t001443 = ('T001443', 0.74*Ms2Mj, 86*pc_in_au, 0.0271, 0.0049, 0, 1e8, 
                825, 2459445, None, None, 
                10.67, 2.2, 'ethraid/data/clean_curves/T001443_2200_AO.csv', [34.4, 0.23*Ms2Mj])
params_219134 = ('219134', 0.79*Ms2Mj, 6.53*pc_in_au, -0.00072, 0.00015, 0, 1e8,
                6349, 2455371, 0.14629, 0.06070, 
                5.56, 0.832, 'ethraid/data/clean_curves/219134_832_speckle.csv')
params_12572 = ('12572', 0.91*Ms2Mj, 65.9*pc_in_au, -0.0608, 0.0042, 5e-8, 6e-6,
                916, 2459320, 0.0748781, 0.0451,
                9.2, 0.832, 'ethraid/data/clean_curves/12572_832_speckle.csv')
params_t001669 = ('T001669', 1.13*Ms2Mj, 111.6*pc_in_au, -0.0226, 0.006, 0, 1e10,
                723, 2459410, None, None, 
                10.22, 0.832, 'ethraid/data/clean_curves/TOI1669_832_speckle.csv')

###########################################################
# T001174 but with idealized constraints from the Vortex coronagraph
params_t001174_vtx = ('T001174_vtx', 0.82*Ms2Mj, 94.9*pc_in_au, -0.116, 0.019, 4.7e-5, 7.2e-5, 
                    748, 2459209, None, None, 
                    10.96, 3.77, 'ethraid/data/clean_curves/vortex_Lband.csv')
             
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
                

# CLS star with known companion.
params_hd186408= ('HD186408', 1065.317938, 4459779.594497, -0.004586, 0.000198, 0, 1e8, 
                12009.069844, 2452962.453122, 0.190556, 0.040519, 840.84, 0.99*Ms2Mj)

# CLS star with known companion from Lea's table 5.
params_hd131156= ('HD131156', 935.993904, 1381916.161366, 0.065991, 0.001904, 0, 1e8, 
                11998.216066, 2452956.829933, 19.315683, 0.175772, 35, 0.65*Ms2Mj)

# CLS star with known companion from Lea's table 5.
params_hd40397= ('HD40397', 925.627738, 4785726.363004, -0.029115, 0.000287, 0, 1e8, 
                8042.117923, 2454859.828823, 2.058092, 0.03366, 60.84, 0.27*Ms2Mj)


################ Marc Hon system HIP73136 aka 8 UMi #################
params_8umi= ('8UMi', 2.26*Ms2Mj, 149.9*pc_in_au, 0.033, 0.018, 0, 1e8, 
                439, 2459570, 0.07669, 0.03175)


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
                
                
                
                
                
                
                