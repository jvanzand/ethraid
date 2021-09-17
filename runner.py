import matplotlib.pyplot as plt
import astropy.constants as c
import numpy as np
# from astropy.time import Time
from scipy.stats import loguniform, beta
import time

import radvel as rv
# import matplotlib.pyplot as plt
# import matplotlib.patches as ptch

from trends import helper_functions_wrapper as hlpw
import plotter
import system_params as sp


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
      

# rv_epoch is the epoch where DATA values of g_dot and g_ddot are computed. Taken from radvel setup file.
m_star, d_star, gammadot, gammadot_err, gammaddot, gammaddot_err,\
        rv_baseline, max_rv, rv_epoch, delta_mu, delta_mu_err = sp.params_synth


# min_per is 4xbaseline because we see ~no curvature yet.
# min_per = 4*rv_baseline
min_per = rv_baseline
min_K = max_rv

min_m = rv.utils.Msini(min_K, min_per, m_star, e=0, Msini_units='jupiter')
min_a = rv.utils.semi_major_axis(min_per, (m_star + min_m*(M_jup/M_sun)))

# # Experimental: revised min_m using gdot. The idea is to find the smallest mass that could produce the observed gdot at the known minimum period. This is not completely right because it uses min_K to get min_m, min_m to get min_a, and then min_a to get a new value for min_m.
#
# min_m = (gammadot*100/(24*3600))*((min_per*24*3600)/6.283185)**2*(m_star*M_sun)/(min_a*14959787070000.0) / M_jup
# print('!!!', min_m)

print('Min m is: ', min_m)
print('Min a is: ', min_a)

# Sampling limits for a and m. Note that if the min_a or min_m parameters fall outside these bounds, the plot will look weird. I can modify later to throw an error, but it's mostly visual.
# # 191939
# # min_a = 0.5
# # min_m = 0.5
# a_lim = (0.8*min_a, 5e1)
# m_lim = (0.8*min_m, 1e2)
# # HIP97166
# a_lim = (1.9, 2e3)
# m_lim = (0.03, 2e3)
# # HD6101
# a_lim = (0.9*min_a, 6e1)
# m_lim = (0.9*min_m, 1e5)
# # HD91204
# a_lim = (0.8*min_a, 8e1)
# m_lim = (0.8*min_m, 1e5)
# synthetic
a_lim = (0.8*min_a, 5e1)
m_lim = (0.8*min_m, 1e2)
# # GL758
# a_lim = (0.5*min_a, 2e2)
# m_lim = (0.5*min_m, 4e2)
# # HIP67246
# a_lim = (0.5*min_a, 2e2)
# m_lim = (0.5*min_m, 4e2)
# # HIP63510
# a_lim = (0.5*min_a, 2e2)
# m_lim = (0.5*min_m, 1e3)
# # 12572
# a_lim = (0.5*min_a, 5e1)
# m_lim = (0.5*min_m, 1e2)
print(a_lim[0], min_a)

grid_num = 100
num_points = int(1e6)
np.set_printoptions(threshold=np.inf)



a_list, m_list, per_list, e_list, i_list, om_list, M_anom_0, E_anom_rv, a_inds, m_inds = \
                                            hlpw.make_arrays(m_star, a_lim, m_lim, rv_epoch, grid_num, num_points)


print('made arrays')

##
start_time = time.time()
##

# Create an array with 1s in allowed regions and 0s in disallowed regions
min_index_m = int(np.ceil(hlpw.value2index(min_m, (0, grid_num-1), m_lim)))
min_index_a = int(np.ceil(hlpw.value2index(min_a, (0, grid_num-1), a_lim)))

prior_array = np.ones((grid_num, grid_num))
prior_array[0:min_index_m, :] = 0
prior_array[:, 0:min_index_a] = 0

# plt.imshow(prior_array, origin='lower')
# plt.show()

# Some targets aren't in the Hip/Gaia catalog, so we can't make the astrometry posterior for them.
try:
    astro_list = hlpw.astro_post(delta_mu, delta_mu_err, m_star, d_star, a_list,
                                 m_list, per_list, e_list, i_list, om_list,
                                 M_anom_0, num_points, grid_num)                      
                                 
    post_astro = hlpw.prob_array(astro_list, a_inds, m_inds, grid_num) * prior_array
    post_astro = post_astro/post_astro.sum()

except:
    astro_list = np.ones(num_points)
    post_astro = np.ones((grid_num, grid_num))
    print('No astrometry data provided. Bounds will be based on RVs only.')
    


rv_list = hlpw.rv_post(gammadot, gammadot_err, gammaddot, gammaddot_err, m_star, 
                        a_list, m_list, per_list, e_list, i_list, om_list, E_anom_rv, 
                        num_points, grid_num)
post_rv = hlpw.prob_array(rv_list, a_inds, m_inds, grid_num) * prior_array
post_tot = hlpw.post_tot(rv_list, astro_list, grid_num, a_inds, m_inds) * prior_array

post_rv = post_rv/post_rv.sum()
post_tot = post_tot/post_tot.sum()

##
end_time = time.time()
##
print('{:.0e} points ran in {:.2f} seconds.'.format(num_points, end_time-start_time))


# plt.imsave('post_rv.png', post_rv, origin='lower')
# plt.imsave('post_astro.png', post_astro, origin='lower')
# plt.imsave('post_tot.png', post_tot, origin='lower')
#
# plt.imshow(post_rv, origin='lower')
# plt.show()
# plt.imshow(post_astro, origin='lower', cmap='jet')
# plt.show()
# plt.imshow(post_tot, origin='lower')
# plt.show()

plotter.joint_plot(m_star, post_tot, post_rv, post_astro, grid_num, a_lim, m_lim, (min_a, min_m), period_lines = True)
print('Plotted')


