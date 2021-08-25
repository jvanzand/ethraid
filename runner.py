import matplotlib.pyplot as plt
import astropy.constants as c
import numpy as np
from astropy.time import Time
from scipy.stats import loguniform, beta

import radvel as rv
import matplotlib.pyplot as plt
import matplotlib.patches as ptch

from trends import helper_functions_wrapper as hlpw


## Constants ##
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30
pc_in_cm = 3.086e18

# params_star = (m_star, distance(cm), gdot, gdot_err, gddot, gddot_err, 
#               rv_baseline(days), max_rv of residuals, rv_epoch, delta_mu, delta_mu_err)

params_191939 = (0.807, 58.3*pc_in_cm, 0.114, 0.006, -6e-5, 1.9e-5, 
                430.2527364352718, 40.0581900021484, 2458847.780463, 0.12767382507786398, 0.034199052901953214)

# HIP97166 params. rv_baseline and max_rv are estimated for ease.
params_HIP97166 = (0.91, 68*pc_in_cm, 0.013, 0.03, -3e-6, 3.2e-5,
                    440, 1, 2458683.353840575, 0.036354497, 0.037699807)

# 12572 params. rv_baseline and max_rv are estimated for ease. gddot := 0, so borrow its error from gdot
params_12572 = (0.91, 65.9*pc_in_cm, -0.0595, 0.0032, 0, 0.0032,
                550, 30, 2458991.236308, 0.0748781, 0.045100458)
                
# HD6106, the highest-trend CLS target with 3sig astro and RV. No curv given. Epoch estimated. Trend/curv taken from manual radvel run.          
params_hd6101 = (0.79, 21.5*pc_in_cm, 0.25858, 0.022, -4.64671e-05, 4.8e-05,
                1626, 300, 2457654.089, 19.411242, 0.240018)
                
# HD91204, a high-trend CLS target with 3sig astro and RV. No curv given. Epoch estimated.           
params_hd91204 = (1.15, 51.8*pc_in_cm, -144.69, 0.914247, 0, 0.1,
                6485.728, 2000, 2455232.024, 3.136, 0.0464) 
                
# T001194 params. rv_baseline and max_rv are estimated for ease.
params_T001194 = (0.98, 150.3*pc_in_cm, 0.019, 0.023, -6.2e-5, 5.8e-5,
                    567, 1, 2458917.385183, None, None)

# Synthetic planet params to test trend code
params_synthetic = (0.79, 21.5*pc_in_cm, 0.059206530869421, 0.00592065308694, -2.0423039552e-4, -2.0423039552e-5,
                1626, 300, 2457654.089, 9.762827453051148, 0.9762827453051148)

# rv_epoch is the epoch where DATA values of g_dot and g_ddot are computed. Taken from radvel setup file.
m_star, d_star, gammadot, gammadot_err, gammaddot, gammaddot_err,\
        rv_baseline, max_rv, rv_epoch, delta_mu, delta_mu_err = params_191939


# min_per is 4xbaseline because we see ~no curvature yet.
min_per = 4*rv_baseline
# min_per = rv_baseline
min_K = max_rv

min_m = rv.utils.Msini(min_K, min_per, m_star, e=0, Msini_units='jupiter')
min_a = rv.utils.semi_major_axis(min_per, (m_star + min_m*(M_jup/M_sun)))

# print('!!!', min_m)
# # Experimental: revised min_m using gdot. The idea is to find the smallest mass that could produce the observed gdot at the known minimum period. This is not completely right because it uses min_K to get min_m, min_m to get min_a, and then min_a to get a new value for min_m.
#
# min_m = (gammadot*100/(24*3600))*((min_per*24*3600)/6.283185)**2*(m_star*M_sun)/(min_a*14959787070000.0) / M_jup
# print('!!!', min_m)

print('Min m is: ', min_m)
print('Min a is: ', min_a)

# Sampling limits for a and m. Note that if the min_a or min_m parameters fall outside these bounds, the plot will look weird. I can modify later to throw an error, but it's mostly visual.
# 191939
min_a = 0.5
min_m = 0.5
a_lim = (0.8*min_a, 5e3)
m_lim = (0.8*min_m, 2e4)
# # HIP97166
# a_lim = (1.9, 2e3)
# m_lim = (0.03, 2e3)
# # HD6101
# a_lim = (0.9*min_a, 6e1)
# m_lim = (0.9*min_m, 1e5)
# # HD91204
# a_lim = (0.8*min_a, 8e1)
# m_lim = (0.8*min_m, 1e5)
# # synthetic
# a_lim = (0.8*min_a, 5e1)
# m_lim = (0.8*min_m, 1e5)
print(a_lim[0], min_a)

grid_num = 100

num_points = int(1e7)

tick_num = 6
tick_size = 30
np.set_printoptions(threshold=np.inf)



a_list, m_list, per_list, e_list, i_list, om_list, M_anom_0, E_anom_rv, a_inds, m_inds = \
                                            hlpw.make_arrays(m_star, a_lim, m_lim, rv_epoch, grid_num, num_points)


print('made arrays')

fig, ax = plt.subplots(figsize=(12,12), dpi = 300)


# Some targets aren't in the Hip/Gaia catalog, so we can't make the astrometry posterior for them.
if delta_mu == None or delta_mu_err == None:
    astro_list = np.ones(num_points)
    post_astro = np.ones((grid_num, grid_num))
    print('all_ones')
    
else:
    astro_list = hlpw.astro_post(delta_mu, delta_mu_err, m_star, d_star, a_list,
                                 m_list, per_list, e_list, i_list, om_list,
                                 M_anom_0, num_points, grid_num)                      
                                 
    post_astro = np.array(hlpw.prob_array(astro_list, a_inds, m_inds, grid_num))
    post_astro = post_astro/post_astro.sum()
    
    t_contours_astro = hlpw.contour_levels(post_astro, [1,2])
    post_astro_cont = ax.contourf(post_astro, t_contours_astro, cmap='Blues', extend='max', alpha=0.5)
    

rv_list = hlpw.rv_post(gammadot, gammadot_err, gammaddot, gammaddot_err, m_star, 
                        a_list, m_list, per_list, e_list, i_list, om_list, E_anom_rv, 
                        num_points, grid_num)
post_rv = np.array(hlpw.prob_array(rv_list, a_inds, m_inds, grid_num))
post_tot = np.array(hlpw.post_tot(rv_list, astro_list, grid_num, a_inds, m_inds))

post_rv = post_rv/post_rv.sum()
post_tot = post_tot/post_tot.sum()

t_contours_rv = hlpw.contour_levels(post_rv, [1,2])
t_contours_tot = hlpw.contour_levels(post_tot, [1,2])

post_rv_cont = ax.contourf(post_rv, t_contours_rv, cmap='Greens', extend='max', alpha=0.5)
post_tot_cont = ax.contourf(post_tot, t_contours_tot, cmap='Reds', extend='max', alpha=0.75)

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


# While the above a_list and m_list are the random samples, these are log-uniform lists for plotting.
a_list = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
m_list = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)


min_index_m = hlpw.value2index(min_m, (0, grid_num-1), m_lim)
min_index_a = hlpw.value2index(min_a, (0, grid_num-1), a_lim)


# Print out the 2-sigma boundaries for the total posterior
bounds = hlpw.bounds_1D(post_tot, [m_lim, a_lim], interp_num = 1e4)
print('a_lim, m_lim = ', bounds[0], bounds[1])


mass_rect = ptch.Rectangle((0, 0), grid_num-1, min_index_m, color='gray', alpha=1.0)
a_rect = ptch.Rectangle((0, 0), min_index_a, grid_num-1, color='gray', alpha=1.0)

ax.add_patch(mass_rect)
ax.add_patch(a_rect)
###################################################

############### Adding labels #####################
label_size = 50
region_label_size = 50
restricted_region_label_size = 40

plt.text((19/32)*grid_num, (7/8)*grid_num, 'RV', size=region_label_size, rotation=50)
plt.text((9/16)*grid_num, (1/4)*grid_num, 'Astrometry', size=region_label_size)

plt.text((1/6)*grid_num, (1/3)*(min_index_m-1), 'Masses disallowed by RVs', size=restricted_region_label_size)
plt.text((1/3)*(min_index_a-1), (1/8)*grid_num, 'Ruled out by minimum period', size=restricted_region_label_size, rotation=90)


ax.set_xlabel('Semi-major Axis (au)', size=label_size)
ax.set_ylabel(r'$M_p$ ($M_{Jup}$)', size=label_size)
###################################################

######## Adding lines of constant period ##########
hip_times  = [Time(1989.85, format='decimalyear').jd, Time(1993.21, format='decimalyear').jd] #https://www.cosmos.esa.int/web/hipparcos/catalogue-summary

gaia_times = [Time('2014-07-25', format='isot').jd, Time('2017-05-28', format='isot').jd] #https://www.cosmos.esa.int/web/gaia/earlydr3

# Time between the midpoints of the two missions
baseline_days = ((gaia_times[1] + gaia_times[0])/2 - (hip_times[1] + hip_times[0])/2)

# Log-spaced masses in Jupiter masses
const_per_m_list = np.logspace(np.log10(min_m), np.log10(m_lim[1]))
const_per_m_inds = hlpw.value2index(const_per_m_list, (0, grid_num-1), m_lim)

for f in range(5):
    
    const_per_a_list = hlpw.period_lines(const_per_m_list, baseline_days/(f+1), m_star)
    const_per_a_inds = hlpw.value2index(const_per_a_list, (0, grid_num-1), a_lim)
    
    values_in_bounds = np.where(const_per_a_list >= min_a)

    plt.plot(const_per_a_inds[values_in_bounds], const_per_m_inds[values_in_bounds], '--k', alpha=0.5)
###################################################

tick_array = np.linspace(0, grid_num-1, tick_num).astype(int)

# Need the .item() method for larger companions because once m_pl > 0.1*M_star, rvutils.Msini() returns [value] instead of value.
plt.xticks(tick_array, [np.round(a_list[i], 1).item() for i in tick_array], size=tick_size)
plt.yticks(tick_array, [np.round(m_list[i], 1).item() for i in tick_array ], size=tick_size)


fig.tight_layout()
fig.savefig('5thCompConstraints_RV_astr.png')
# fig.savefig('plots/5thCompConstraints_RV_astr.png')


##########################################
plt.close()
fig, ax = plt.subplots(1,2, figsize=(12,12))
sma_1d = post_tot.sum(axis=0)
mass_1d = post_tot.sum(axis=1)

ax[0].plot(range(grid_num), sma_1d)
plt.sca(ax[0])
plt.xticks(tick_array, [np.round(a_list[i], 1).item() for i in tick_array], size=10)

ax[1].plot(range(grid_num), mass_1d)
plt.sca(ax[1])
plt.xticks(tick_array, [np.round(m_list[i], 1).item() for i in tick_array], size=10)

fig.savefig('1_d_posts.png')
# fig.savefig('plots/1_d_posts.png')

############################################
