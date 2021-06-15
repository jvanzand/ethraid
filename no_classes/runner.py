import matplotlib.pyplot as plt
import astropy.constants as c
import numpy as np

# replace with hlpw
# import helper_functions_wrapper as hlp


import radvel as rv

from c_kepler import _kepler as ck

from scipy.stats import loguniform, beta

# import post_functions_python as  pfp
import post_functions as  pfp

## Constants ##

pi = 3.141592653589793
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30


gammadot      = 0.114
gammadot_err  = 0.006
gammaddot     = -6e-5
gammaddot_err = 1.9e-5

delta_mu     = 0.12767382507786398
delta_mu_err = 0.034199052901953214


m_star = 0.807
a_lim = (1.9, 5e1)
m_lim = (1.5, 2e2)
grid_num = 30
num_points = int(1e4) 
t_num = 100
tick_num = 6
tick_size = 30


a_list, m_list, per_list, e_list, i_list, om_list, M_anom_list, E_anom_list, T_anom_list, a_inds, m_inds = pfp.make_arrays(m_star, a_lim, m_lim, grid_num, num_points)
                                                
print('made arrays')

post_rv    = pfp.rv_post(gammadot, gammadot_err, gammaddot, gammaddot_err, m_star, a_list, m_list, per_list, e_list, i_list, om_list, E_anom_list, num_points, grid_num, a_inds, m_inds)

# post_astro =                        pfp.astro_post(delta_mu, delta_mu_err, m_star, a_list, m_list, per_list, e_list, i_list, om_list, T_anom_list, num_points, grid_num, a_inds, m_inds, t_num)














# import matplotlib.pyplot as plt
# importmatplotlib.patches as ptch

# plt.imshow(post_astro, origin='lower')
# plt.show()
# # The priors for minimum period and planet mass. min_per is 4xbaseline because we see ~no curvature yet.
# rv_baseline = 430.2527364352718
# max_rv = 40.0581900021484
# min_per = 4*rv_baseline
#
# # While the above a_list and m_list are the random samples, these are log-uniform lists for plotting.
# a_list = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
# m_list = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)
#
# min_m = rv.utils.Msini(max_rv, min_per, m_star, e=0, Msini_units='jupiter')
# min_a = rv.utils.semi_major_axis(min_per, (m_star + min_m*(M_jup/M_sun)))
#
#
# min_index_m = hlp.value2index(min_m, (0, grid_num-1), m_lim)
# min_index_a = hlp.value2index(min_a, (0, grid_num-1), a_lim)
#
#
# t_contours_rv = hlp.contour_levels(post_rv, [1,2])
#
# bounds = hlp.bounds_1D(post_rv, [m_lim, a_lim], interp_num = 1e4)
# print('a_lim, m_lim = ', bounds[0], bounds[1])
#
# """
# fig, ax = plt.subplots(figsize=(12,12))
#
# post_rv_cont = ax.contourf(post_rv, t_contours_rv, cmap='Greens', extend='max', alpha=0.5)
#
# mass_rect = ptch.Rectangle((0, 0), grid_num-1, min_index_m, color='gray', alpha=1.0)
# a_rect = ptch.Rectangle((0, 0), min_index_a, grid_num-1, color='gray', alpha=1.0)
#
# ax.add_patch(mass_rect)
# ax.add_patch(a_rect)
#
# tick_array = np.linspace(0, grid_num-1, tick_num).astype(int)
#
# plt.xticks(tick_array, [np.round(a_list[i], 1) for i in tick_array], size=tick_size)
# plt.yticks(tick_array, [np.round(m_list[i], 1) for i in tick_array ], size=tick_size)
#
#
# fig.tight_layout()
# # fig.savefig('5thCompConstraints_RV_astr.png')
# plt.show()
# """