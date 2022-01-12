M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
Ms2Mj = M_sun/M_jup

# params_star = (star_name, m_star, distance(cm), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), rv_epoch, delta_mu, delta_mu_err)

############## Syetems for papers #########################
params_191939 = ('191939', 0.807*Ms2Mj, 58.3*pc_in_au, 0.1116, 0.0037, -3.44e-5, 5.1e-6, 
                778.855, 2459192.641, 0.12767382507786398, 0.034199052901953214)

# Potential interesting DG paper. No gddot defined in RVS preferred fit.
params_12572 = ('12572', 0.91*Ms2Mj, 65.9*pc_in_au, -0.0613, 0.0025, 0, 0.1,
                600, 2459183.889, 0.075, 0.045)
###########################################################

# CLS system in automation. Testing locally to debug.
params_hd24916 = ('HD24916_local', 757.9, 3252875.04, -0.0038566, 0.000666, 0, 1e10,
                    8110.823, 2454421.49, 0.4271376, 0.052581)

params_191939_old = ('191939_old', 0.807*Ms2Mj, 58.3*pc_in_au, 0.114, 0.006, -6e-5, 1.9e-5, 
                430.2527364352718, 2458847.780463, 0.12767382507786398, 0.034199052901953214)


# HIP97166 params. rv_baseline and max_rv are estimated for ease.
params_HIP97166 = (0.91*Ms2Mj, 68*pc_in_au, 0.013, 0.03, -3e-6, 3.2e-5,
                    440, 2458683.353840575, 0.036354497, 0.037699807)
                
# HD6106, the highest-trend CLS target with 3sig astro and RV. No curv given. Epoch estimated. Trend/curv taken from manual radvel run.          
params_hd6101 = (0.79*Ms2Mj, 21.5*pc_in_au, 0.25858, 0.022, -4.64671e-05, 4.8e-05,
                1626, 2457654.089, 19.411242, 0.240018)
                
# HD91204, a high-trend CLS target with 3sig astro and RV. No curv given. Epoch estimated.           
params_hd91204 = (1.15*Ms2Mj, 51.8*pc_in_au, -144.69, 0.914247, 0, 0.1,
                6485.728, 2455232.024, 3.136, 0.0464)
                
# HD238894, Paul Dalba's target with no HIP ID. Can only be constrained using RVs.
# gammas taken from Table 3 in paper. good agreement with preferred radvel       
params_hd238894 = (1.16*Ms2Mj, 117.6*pc_in_au, -0.1205, 0.0043, 2.14e-4, 3.9e-5,
                570, 2459177.6, None, None)




################ Validation Systems #################
         
# GL758, an example star in Tim Brandt's Orvara code. Using this to compare results.
params_gl758 = (0.95*Ms2Mj, 15.5*pc_in_au, -0.00633, 0.00025, -8.19e-7, 0.67e-7,
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
params_synth = ('synthetic', 1249.7, 15267565.22, -0.861595, 0.005, -8.58e-6, 10000000,
                7508.3, 2454737.94, 0.16233, 0.0965)

# # # Synthetic: a = 17, m = 47, i = pi/2, tiny astro signal
# params_synth = (1*Ms2Mj, 10*pc_in_au, -0.02944, 0.02944/5, 1.88e-5, 1.88e-5/5,
#                 300, 2459000, 0.59, 0.59/5)
                
# # # Synthetic: a = 17, m = 47, i = pi/3, larger astro signal
# params_synth = (1*Ms2Mj, 10*pc_in_au, -0.0255, 0.0231/5, 1.63e-5, 1.63e-5/5,
#                 300, 2459000, 3.282, 3.283/5)
                
                
                
                
                
                
                