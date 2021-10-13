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
                
# HD6106, Legacy detection. Binary with P ~ 10k day/29 year, M_comp = 1.17-0.79 = 0.38 (Balega2006)   
params_hd6101 = (0.79, 21.5*pc_in_cm, 0.254547, 0.0211, -1.008e-4, 8.825e-05,
                1626, 1000, 2457663.942, 19.411242, 0.240018)
                
# HD91204, a high-trend CLS target with 3sig astro and RV. No curv given. Epoch estimated.           
params_hd91204 = (1.15, 51.8*pc_in_cm, -144.69, 0.914247, 0, 0.1,
                6485.728, 2000, 2455232.024, 3.136, 0.0464)
                
# HD238894, Paul Dalba's target with no HIP ID. Can only be constrained using RVs.
# gammas taken from Table 3 in paper. good agreement with preferred radvel       
params_hd238894 = (1.16, 117.6*pc_in_cm, -0.1205, 0.0043, 2.14e-4, 3.9e-5,
                570, 40, 2459177.6, None, None)




################ Validation Systems #################
         
# GL758, an example star in Tim Brandt's Orvara code. Using this to compare results.
params_gl758 = (0.95, 15.5*pc_in_cm, -0.00633, 0.00025, -8.19e-7, 0.67e-7,
                8413.010, 60, 2454995.123, 1.0397, 0.0261)
                
# HD4747, an example star in Tim Brandt's Orvara code. Using this to compare results.
params_hd4747 = (0.83, 18.8*pc_in_cm, -0.10605, 0.10605/5, 1.1276e-5, 1.1276e-5/5,
                8473, 200, 2455015.622, 3.74, 0.049)
                
# HIP67246, a confirmed giant published by Sarah Blunt. I truncated the timeseries to get gdot/gddot as a test.
# gdot/gddot fit with scipy.optimize.curve_fit bc rvsearch gave weird results
params_hip67246 = (1.10, 30.6*pc_in_cm, 0.01, 0.002, 4.53e-5, 0.72e-5,
                1126.94, 50, 2457685.312, 0.135, 0.035)
         
# HIP63510 aka Ross 458, a binary star system (primary is 0.6 M_sun, secondary is 80 M_J.)
# https://arxiv.org/pdf/1103.3544.pdf
params_hip63510 = (0.6, 11.4*pc_in_cm, 1.53, 3*0.153, -1.9e-4, 1.9e-5,
                2336.7, 2, 2455715.384, 14.465, 0.23)
                

# # HD 3795, Legacy detection. Jump says low-mass star ~300 M_J at ~20 AU
# params_hd3795 = (0.85, 28.6*pc_in_cm, 0.327, 0.0019, -3.226e-5, 2.63e-6,
#                 5135.03, 2000, 2452933.471, 23.8876, 0.089482)
#
# # HD 3795 again, but now with full curv (not just trend)
# params_hd3795_curv = (0.85, 28.6*pc_in_cm, 0.6585, 0.01155, 7.082e-5, 4.283e-6,
#                 9118, 2000, 2454924.91172, 23.8876, 0.089482)

# HD 4614, Legacy trend detection.
params_hd4614 = (0.91, 6.0*pc_in_cm, 0.004466, 0.001344, -2.153e-6, 4.79e-7,
                9505, 100, 2454722.410, 6.01, 0.19523)

# HD 9446, Legacy warm Jupiter at 190 days (second giant at 30 days). Truncated for trend.
params_9446 = (1.05, 52.9*pc_in_cm, 1.18499, 0.21784, 0.00407164, 0.02181313,
                88.9, 45, 2458758.475, 0.10379, 0.05413)
                
# HD 9986, Legacy target showing ~3000 day planet. Jump comments suggest possible activity. Truncated.
params_9986 = (1.03, 25.7*pc_in_cm, -0.0275133, 0.00315776, 4.30867e-5, 1.8516651e-5,
                1307.3, 15, 2458803.4015, 0.03779, 0.052464)


#######Synthetic Experiment: fix m = 30, a = 5; vary e, i, om, M_anom_0########
###############################################################################
# # Synthetic: m = 30, a = 5, M_0 = 0, om = 0, i = pi/2, e = 0
# params_synth = (1, 10*pc_in_cm, 0.54945, 0.54945/5, 0.00038, 0.00038/3,
#                 300, 250, 2459000, 4.52623, 4.52623/10)

# # Synthetic: m = 30, a = 5, M_0 = 3, om = 2, i = pi/2, e = 0.1
# params_synth = (1, 10*pc_in_cm, 0.30567, 0.30567/5, 0.00111, 0.00111/3,
#                 300, 250, 2459000, 5.95431, 5.95431/10)

# # Synthetic: m = 30, a = 5, M_0 = 2.5*pi, om = 1, i = 2/3 * pi, e = 0.4
# params_synth = (1, 10*pc_in_cm, -1.27110, 1.27110/5, 0.00070, 0.00070/3,
#                 300, 250, 2459000, 4.91312, 4.91312/10)

# # Synthetic: m = 30, a = 5, M_0 = 2.5*pi, om = 1, i = 2/3 * pi, e = 0.0
# params_synth = (1, 10*pc_in_cm, -0.46783, 0.46783/5, -0.00014, 0.00014/3,
#                 300, 250, 2459000, 6.47812, 6.47812/10)

# # Synthetic: m = 30, a = 5, M_0 = 2.5*pi, om = 3pi/5, i = pi/4, e = 0.25
# params_synth = (1, 10*pc_in_cm, -0.61735, 0.61735/5, 0.00143, 0.00143/3,
#                 300, 250, 2459000, 4.30799, 4.30799/10)

# # Synthetic: m = 30, a = 5, M_0 = 0, om = 3pi/5, i = pi/4, e = 0.7
# params_synth = (1, 10*pc_in_cm, 0.15440, 0.15440/5, 0.00024, 0.00024/3,
#                 300,250, 2459000, 3.33498, 3.33498/10)
###############################################################################
######## Sub-experiment: fix m = 15, a = 3.5, om = 3pi/5, i = pi/4, e = 0.1. Vary only M_0

# # M_0 = 0
# params_synth = (1, 10*pc_in_cm, 0.05212, 0.05212/5, -0.00068, 0.00068/3,
#               300, 50, 2459000, 3.37678, 3.37678/10)

# # M_0 = 40 deg
# params_synth = (1, 10*pc_in_cm, -0.23265, 0.23265/5, -0.00091, 0.00091/3,
#               300, 50, 2459000, 2.59227, 2.59227/10)

# # M_0 = 80 deg
# params_synth = (1, 10*pc_in_cm, -0.49059, 0.49059/5, -0.00015, 0.00015/3,
#               300, 50, 2459000, 2.78150, 2.78150/10)

# # M_0 = 120 deg
# params_synth = (1, 10*pc_in_cm, -0.43963, 0.43963/5, 0.00134, 0.00134/3,
#               300, 50, 2459000, 3.33861, 3.33861/10)

# # M_0 = 160 deg
# params_synth = (1, 10*pc_in_cm, -0.10003, 0.10003/5, 0.00162, 0.00162/3,
#               300, 50, 2459000, 3.40860, 3.40860/10)

# # M_0 = 200 deg
# params_synth = (1, 10*pc_in_cm, 0.19506, 0.19506/5, 0.00099, 0.00099/3,
#               300, 50, 2459000, 3.01655, 3.01655/10)

# # M_0 = 240 deg
# params_synth = (1, 10*pc_in_cm, 0.32331, 0.32331/5, 0.00044, 0.00044/3,
#               300, 50, 2459000, 2.52898, 2.52898/10)

# # M_0 = 280 deg
# params_synth = (1, 10*pc_in_cm, 0.32628, 0.32628/5, 0.00007, 0.00007/3,
#               300, 50, 2459000, 2.71057, 2.71057/10)

# # M_0 = 320 deg
# params_synth = (1, 10*pc_in_cm, 0.23752, 0.23752/5, -0.00027, 0.00027/3,
#               300, 50, 2459000, 3.44293, 3.44293/10)
              
              
              
              
              
              
                