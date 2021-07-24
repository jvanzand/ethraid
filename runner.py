import matplotlib.pyplot as plt
import astropy.constants as c
import numpy as np
from scipy.stats import loguniform, beta

import radvel as rv

from trends import helper_functions_wrapper as hlpw


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
pc_in_cm = 3.086e18

params_191939 = (0.807, 1.7989500299953727e+20, 0.114, 0.006, -6e-5, 1.9e-5, 
                430.2527364352718, 40.0581900021484, 2458847.780463, 0.12767382507786398, 0.034199052901953214)

# HIP97166 params. rv_baseline and max_rv are estimated for ease.
params_HIP97166 = (0.91, 68*pc_in_cm, 0.013, 0.03, -3e-6, 3.2e-5,
                    440, 1, 2458683.353840575, 0.036354497, 0.037699807)

# 12572 params. rv_baseline and max_rv are estimated for ease.
params_12572 = (0.91, 65.9*pc_in_cm, -0.0595, 0.0032, 0, 0.0032,
                550, 30, 2458991.236308, 0.0748781, 0.045100458)

# rv_epoch is the epoch where DATA values of g_dot and g_ddot are computed. Taken from radvel setup file.
m_star, d_star, gammadot, gammadot_err, gammaddot, gammaddot_err,\
        rv_baseline, max_rv, rv_epoch, delta_mu, delta_mu_err = params_HIP97166

# Sampling limits for a and m. Note that if the min_a or min_m parameters fall outside these bounds, the plot will look weird. I can modify later to throw an error, but it's mostly visual.
# a_lim = (1.9, 5e1)
# m_lim = (1.5, 2e2)
a_lim = (1.9, 2e2)
m_lim = (0.03, 2e2)

grid_num = 100

num_points = int(1e4)

t_num = 2
tick_num = 6
tick_size = 30
np.set_printoptions(threshold=np.inf)


a_list, m_list, per_list, e_list, i_list, om_list, E_anom_rv, T_anom_astro, a_inds, m_inds = \
                                            hlpw.make_arrays(m_star, a_lim, m_lim, rv_epoch, grid_num, num_points)


print('made arrays')

rv_list = hlpw.rv_post(gammadot, gammadot_err, gammaddot, gammaddot_err, m_star, 
                        a_list, m_list, per_list, e_list, i_list, om_list, E_anom_rv, 
                        num_points, grid_num)
astro_list = hlpw.astro_post(delta_mu, delta_mu_err, m_star, d_star, a_list,
                             m_list, per_list, e_list, i_list, om_list,
                             T_anom_astro, num_points, grid_num, t_num)

# np.save('erik_samples/rv_probs.npy', rv_list)
# np.save('erik_samples/astro_probs.npy', astro_list)
# np.save('erik_samples/a_indices.npy', a_inds)
# np.save('erik_samples/m_indices.npy', m_inds)
#
# rv_list = np.load('erik_samples/rv_probs.npy')
# astro_list = np.load('erik_samples/astro_probs.npy')
                                                               
post_rv = np.array(hlpw.prob_array(rv_list, a_inds, m_inds, grid_num))
post_astro = np.array(hlpw.prob_array(astro_list, a_inds, m_inds, grid_num))
post_tot = np.array(hlpw.post_tot(rv_list, astro_list, grid_num, a_inds, m_inds))


post_rv = post_rv/post_rv.sum()
post_astro = post_astro/post_astro.sum()
post_tot = post_tot/post_tot.sum()



import matplotlib.pyplot as plt
import matplotlib.patches as ptch
#
# plt.imsave('post_rv.png', post_rv, origin='lower')
# plt.imsave('post_astro.png', post_astro, origin='lower')
# plt.imsave('post_tot.png', post_tot, origin='lower')
#
# plt.imshow(post_rv, origin='lower')
# plt.show()
# plt.imshow(post_astro, origin='lower')
# plt.show()
# plt.imshow(post_tot, origin='lower')
# plt.show()


# min_per is 4xbaseline because we see ~no curvature yet.
min_per = 4*rv_baseline

# While the above a_list and m_list are the random samples, these are log-uniform lists for plotting.
a_list = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
m_list = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)

min_m = rv.utils.Msini(max_rv, min_per, m_star, e=0, Msini_units='jupiter')
min_a = rv.utils.semi_major_axis(min_per, (m_star + min_m*(M_jup/M_sun)))


min_index_m = hlpw.value2index(min_m, (0, grid_num-1), m_lim)
min_index_a = hlpw.value2index(min_a, (0, grid_num-1), a_lim)


t_contours_rv = hlpw.contour_levels(post_rv, [1,2])
t_contours_astro = hlpw.contour_levels(post_astro, [1,2])
t_contours_tot = hlpw.contour_levels(post_tot, [1,2])

# bounds = hlpw.bounds_1D(post_tot, [m_lim, a_lim], interp_num = 1e4)
# print('a_lim, m_lim = ', bounds[0], bounds[1])


fig, ax = plt.subplots(figsize=(12,12), dpi = 300)

post_astro_cont = ax.contourf(post_astro, t_contours_astro, cmap='Blues', extend='max', alpha=0.5)
post_rv_cont = ax.contourf(post_rv, t_contours_rv, cmap='Greens', extend='max', alpha=0.5)
post_tot_cont = ax.contourf(post_tot, t_contours_tot, cmap='Reds', extend='max', alpha=0.75)

mass_rect = ptch.Rectangle((0, 0), grid_num-1, min_index_m, color='gray', alpha=1.0)
a_rect = ptch.Rectangle((0, 0), min_index_a, grid_num-1, color='gray', alpha=1.0)

ax.add_patch(mass_rect)
ax.add_patch(a_rect)

###################################################
label_size = 50
region_label_size = 50
restricted_region_label_size = 40

plt.text((19/32)*grid_num, (7/8)*grid_num, 'RV', size=region_label_size, rotation=50)
plt.text((9/16)*grid_num, (1/4)*grid_num, 'Astrometry', size=region_label_size)

plt.text((1/6)*grid_num, (1/3)*(min_index_m-1), 'Masses disallowed by RVs', size=restricted_region_label_size)
plt.text((1/3)*(min_index_a-1), (1/8)*grid_num, 'Ruled out by minimum period', size=restricted_region_label_size, rotation=90)


ax.set_xlabel('Semi-major Axis (au)', size=label_size)
ax.set_ylabel(r'$M_p$ ($M_{Jup}$)', size=label_size)
# ax.set_title('Astrometry', size=title_size)
###################################################

tick_array = np.linspace(0, grid_num-1, tick_num).astype(int)

plt.xticks(tick_array, [np.round(a_list[i], 1) for i in tick_array], size=tick_size)
plt.yticks(tick_array, [np.round(m_list[i], 1) for i in tick_array ], size=tick_size)


fig.tight_layout()
fig.savefig('5thCompConstraints_RV_astr.png')
# plt.show()
