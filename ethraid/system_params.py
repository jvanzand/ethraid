M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
Ms2Mj = M_sun/M_jup

# params_star = (star_name, m_star, distance(AU), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), rv_epoch, delta_mu, delta_mu_err, 
#               vmag=None, imag_wavelength=None, contrast_str=None, scatter_tuple=[sma, mass])


# # From Lubin et al. 2022 (params from ~mid-2021)
params_191939 = ('191939', 0.807*Ms2Mj, 53.8*pc_in_au, 0.114, 0.006, -6e-5, 1.9e-5,
                430.2527364352718, 2458847.780463, 0.12767382507786398, 0.034199052901953214,
                8.97)