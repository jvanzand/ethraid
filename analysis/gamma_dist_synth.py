"""
Our main goal has been to be able to take in a set of (gdot, gdotdot) pairs and use this to model a distribution of long-period planets.
This module is trying to do the reverse: assume a distribution of planets and find the resulting (gdot, gdotdot) values that they
would produce. By tweaking the planet distribution, we can try to match a known set of (gdot, gdotdot) pairs and guess the planet
distribution.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FormatStrFormatter

from mpl_toolkits.axes_grid1 import make_axes_locatable
import helper_functions as hlp
import astropy.constants as c
import radvel as rv

from scipy.stats import loguniform, beta
from constants import *

np.random.seed(0)

def value2index_lin(value, index_space, value_space):
    """
    Take a value on a linear scale and convert it to an index. index_space and value_space are expected
    as tuples of the form (min_value, max_value).
    """
    
    value = np.array(value)

    
    index_range = index_space[1] - index_space[0]
    value_range = value_space[1] - value_space[0]
    
    index = (value - value_space[0]) * (index_range/value_range) + index_space[0]
    
    if isinstance(index, float):
        return int(index)
    else:
        return [int(np.round(i)) for i in index]
def gamma_array(a_lim, m_lim, num_points = int(1e3),  grid_num = 200, precision = 1, ecc = False):

    m_star = 1
    a_min, a_max = a_lim
    m_min, m_max = m_lim
    
    # a_list = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
    # m_list = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)
    #
    # gamma_dot_array = np.zeros((len(m_list), len(a_list)))
    # gamma_dotdot_array = np.zeros((len(m_list), len(a_list)))
    
    ######################
    a_list = loguniform.rvs(a_min, a_max, size=num_points)
    m_list= loguniform.rvs(m_min, m_max, size=num_points)
    per_list = hlp.P(a_list, (m_star+m_list*(c.M_jup/c.M_sun).value))
    
    if ecc == True:
        e_list = beta(0.867, 3.03).ppf(np.random.uniform(0, 0.99, num_points))
    elif ecc == False:  
        e_list = np.linspace(0, 0, num_points)
        
    M_anom_list = np.random.uniform(0, 2*np.pi, num_points)
    
    cosi_list = np.random.uniform(0, 1, num_points)
    i_list = np.arccos(cosi_list)
    om_list = np.random.uniform(0, 2*np.pi, num_points)
    
    
    gamma_dot_vec, gamma_dotdot_vec = hlp.gamma(a_list, m_list, per_list, e_list, i_list, om_list, M_anom_list)
    gamma_dot_vec = gamma_dot_vec*365
    gamma_dotdot_vec = gamma_dotdot_vec*365*365
    
    # There are lots of outliers, but we need to set symmetric limits for plotting and consistency
    # I'm choosing +/- 1 standard deviation from the mean of both gdot and gddot
    # fac_gd = 1e-3
    # fac_gdd = 1e-6
    fac_gd = 0.3
    fac_gdd = 0.0025
    global gd_low, gd_high
    global gdd_low, gdd_high
    print(np.median(gamma_dot_vec), np.median(gamma_dotdot_vec))
    # gd_low, gd_high = gamma_dot_vec.mean() - fac*gamma_dot_vec.std(), gamma_dot_vec.mean() + fac*gamma_dot_vec.std()
    # gdd_low, gdd_high = gamma_dotdot_vec.mean() - fac*gamma_dotdot_vec.std(), gamma_dotdot_vec.mean() + fac*gamma_dotdot_vec.std()
    gd_low, gd_high = 0 - fac_gd*gamma_dot_vec.std(), 0 + fac_gd*gamma_dot_vec.std()
    gdd_low, gdd_high = 0 - fac_gdd*gamma_dotdot_vec.std(), 0 + fac_gdd*gamma_dotdot_vec.std()
    
    ########################
    global gd_mean, gdd_mean
    gd_mean, gdd_mean = gamma_dot_vec.mean(), gamma_dotdot_vec.mean()
    ########################
    
    
    # These are the indices in each list where gdot/gddot was outside the 1-sigma range
    gd_cut_inds = list(np.where(gamma_dot_vec < gd_low)[0]) + list(np.where(gd_high < gamma_dot_vec)[0])
    gdd_cut_inds = list(np.where(gamma_dotdot_vec < gdd_low)[0]) + list(np.where(gdd_high < gamma_dotdot_vec)[0])
    
    # I ADD the index lists together. This just means that a point will get dropped if EITHER gdot or gddot is out of range
    gd = np.delete(gamma_dot_vec, gd_cut_inds+gdd_cut_inds)
    gdd = np.delete(gamma_dotdot_vec, gd_cut_inds+gdd_cut_inds)
    
    
    # Linear bins on the interval (mean +/- st_dev) for both gdot and gddot
    global gd_bins, gdd_bins
    gd_bins = np.linspace(gd_low, gd_high, int(grid_num))
    gdd_bins = np.linspace(gdd_low, gdd_high, int(grid_num))
    
    
    # Assign each point an index in the x-direction and the y-direction.
    # digitize has a quirk of bin indexing, so the commented lines below could help if needed
    gd_inds = np.digitize(gd, bins = gd_bins)
    # gd_inds = np.where(gd_inds == grid_num, gd_inds - 1, gd_inds)
    

    gdd_inds = np.digitize(gdd, bins = gdd_bins)
    # gdd_inds = np.where(gdd_inds == grid_num, gdd_inds - 1, gdd_inds)
    
    # New 2D array to bin up the gdot and gddot values
    gamma_array = np.zeros((grid_num, grid_num))
    for i in range(len(gd)):
        gd_i = gd_inds[i]
        gdd_i = gdd_inds[i]
        gamma_array[gdd_i, gd_i] += 1


    return gamma_array
    # return gd, gdd
    # return gamma_dot_vec, gamma_dotdot_vec

tick_size = 30
label_size = 30
title_size = 35
tick_num = 4
#
num_points = int(1e5)
grid_num = int(500)
a_lim = (1.9, 1e2)
m_lim = (1.5, 1e2)


##############
fig, axs = plt.subplots(figsize=(12,10))
fig.tight_layout(pad=10)
# data_no_ecc = gamma_array(a_lim, m_lim, num_points, grid_num, 1, ecc=False)
data_ecc = gamma_array(a_lim, m_lim, num_points, grid_num, 1, ecc=True)

data_ecc = data_ecc/data_ecc.sum()

trend_df = pd.read_csv('../data/trend_list.csv')
trend_df_nonan = trend_df.dropna(subset = ['dvdt', 'curv'])[['hostname', 'dvdt', 'curv']]

# Cast the values into indices for plotting
trend_df_nonan['dvdt_ind'] = value2index_lin(trend_df_nonan['dvdt']*365, (0, grid_num-1), (gd_low, gd_high)) # Convert from m/s/d to m/s/yr
trend_df_nonan['curv_ind'] = value2index_lin(trend_df_nonan['curv']*365**2, (0, grid_num-1), (gdd_low, gdd_high))



tick_array = np.linspace(0, grid_num-1, tick_num).astype(int)



im = axs.imshow(data_ecc, origin='lower')#, cmap='jet')#, norm=LogNorm())
axs.scatter(trend_df_nonan.dvdt_ind, trend_df_nonan.curv_ind, c='red', s=10)


axs.set_xlabel(r'$\dot{{\gamma}}$ (m/s/yr)', size = label_size)
axs.set_ylabel(r'$\ddot{{\gamma}} (m/s/yr^2)$', size = label_size)
axs.set_title(r'$\gamma$ Distribution', size = title_size)

axs.set_xticks(tick_array)
axs.set_xticklabels([np.round(gd_bins[i], 2) for i in tick_array], size = tick_size)
axs.set_yticks(tick_array)
axs.set_yticklabels(['{:.1f}'.format(np.round(gdd_bins[i], 10)) for i in tick_array], size = tick_size)
cb = fig.colorbar(im)
cb.ax.tick_params(labelsize=25)


fig.savefig('gamma_dist.png')
plt.show()















