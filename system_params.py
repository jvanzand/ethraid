pc_in_cm = 3.086e18

# params_star = (m_star, distance(cm), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), max_rv of residuals, rv_epoch, delta_mu, delta_mu_err)


############## Syetems for papers #########################
params_191939 = (0.807, 58.3*pc_in_cm, 0.114, 0.006, -6e-5, 1.9e-5, 
                430.2527364352718, 40.0581900021484, 2458847.780463, 0.12767382507786398, 0.034199052901953214)

# Potential interesting DG paper. No gddot defined in RVS preferred fit.
params_12572 = (0.91, 65.9*pc_in_cm, -0.0613, 0.0025, 0, 0.1,
                600, 20, 2459183.889, 0.075, 0.045)

# params_12572 = (0.91, 65.9*pc_in_cm, -0.1613, 0.0025, 0, 0.1,
#                 600, 10, 2459183.889, 0.075, 0.000045)
###########################################################


# HIP97166 params. rv_baseline and max_rv are estimated for ease.
params_HIP97166 = (0.91, 68*pc_in_cm, 0.013, 0.03, -3e-6, 3.2e-5,
                    440, 1, 2458683.353840575, 0.036354497, 0.037699807)
                
# HD6106, the highest-trend CLS target with 3sig astro and RV. No curv given. Epoch estimated. Trend/curv taken from manual radvel run.          
params_hd6101 = (0.79, 21.5*pc_in_cm, 0.25858, 0.022, -4.64671e-05, 4.8e-05,
                1626, 300, 2457654.089, 19.411242, 0.240018)
                
# HD91204, a high-trend CLS target with 3sig astro and RV. No curv given. Epoch estimated.           
params_hd91204 = (1.15, 51.8*pc_in_cm, -144.69, 0.914247, 0, 0.1,
                6485.728, 2000, 2455232.024, 3.136, 0.0464)
                
# HD238894, Paul Dalba's target with no HIP ID. Can only be constrained using RVs.
# gdot taken from 2planets_inner_transit setup file (no gddot defined).         
params_hd238894 = (1.16, 117.6*pc_in_cm, -0.138, 0.013, 0, 0.1,
                570, 40, 2459177.6, None, None)

# Synthetic planet params to test trend code
params_synthetic = (0.79, 21.5*pc_in_cm, 0.088372596, 0.0088372596, -0.0001712052, -0.00001712052,
                1626, 300, 2457654.089, 47.176, 4.7176)




################ Validation Systems #################
         
# GL758, an example star in Tim Brandt's Orvara code. Using this to compare results.
params_gl758 = (0.95, 15.5*pc_in_cm, -0.00633, 0.00025, -8.19e-7, 0.67e-7,
                8413.010, 60, 2454995.123, 1.0397, 0.0261)
                
# HIP67246, a confirmed giant published by Sarah Blunt. I truncated the timeseries to get gdot/gddot as a test.
# gdot/gddot fit with scipy.optimize.curve_fit bc rvsearch gave weird results
params_hip67246 = (1.10, 30.6*pc_in_cm, 0.01, 0.002, 4.53e-5, 0.72e-5,
                1126.94, 50, 2457685.312, 0.135, 0.035)
         
# HIP63510 aka Ross 458, a binary star system (primary is 0.6 M_sun, secondary is 80 M_J.)
# https://arxiv.org/pdf/1103.3544.pdf
params_hip63510 = (0.6, 11.4*pc_in_cm, 1.53, 3*0.153, -1.9e-4, 1.9e-5,
                2336.7, 500, 2455715.384, 14.465, 0.23)
                
                
                
                
                
                