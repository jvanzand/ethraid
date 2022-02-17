M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
Ms2Mj = M_sun/M_jup

# params_star = (star_name, m_star, distance(au), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), rv_epoch, delta_mu, delta_mu_err)

############## Syetems for papers #########################
params_191939 = ('191939', 0.807*Ms2Mj, 58.3*pc_in_au, 0.1116, 0.0037, -3.44e-5, 5.1e-6, 
                778.855, 2459192.641, 0.12767382507786398, 0.034199052901953214)

# Potential interesting DG paper.
params_12572 = ('12572', 949.1484713, 13587931.9002, -0.0575, 0.0053, -4e-6, 1e-5,
                697.8749, 2459218.7643, 0.0748781, 0.0451)
###########################################################

params_191939_old = ('191939_old', 0.807*Ms2Mj, 58.3*pc_in_au, 0.114, 0.006, -6e-5, 1.9e-5, 
                430.2527364352718, 2458847.780463, 0.12767382507786398, 0.034199052901953214)

# CLS star with known companion from Lea's table 5.
params_hd182488= ('HD182488', 0.96*Ms2Mj, 3.196e6, -0.005390, 0.000345, -9.5e-7, 9e-8, 
                7951.02, 2454980e6, 1.0397, 0.0261)
                
# CLS star with known companion from Lea's table 5.
params_hd201091= ('HD201091', 0.64*Ms2Mj, 7.182422e5, -0.007369, 0.000355, -3.22e-7, 5.4e-8, 
                11886.71, 2452901e6, 4.763910, 0.240117)

# CLS star with known companion from Lea's table 5.
params_hd131156= ('HD131156', 935.993904, 1381916.161366, 0.065991, 0.001904, 0, 1e8, 
                11998.216066, 2452956.829933, 19.315683, 0.175772)

# CLS star with known companion from Lea's table 5.
params_hd40397= ('HD40397', 925.627738, 4785726.363004, -0.029115, 0.000287, 0, 1e8, 
                8042.117923, 2454859.828823, 2.058092, 0.03366)

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
                
                
                
                
                
                
                