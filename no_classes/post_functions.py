from kern_profiler_dummy import *

import matplotlib.pyplot as plt
import astropy.constants as c
import numpy as np


import sys
sys.path.append('../')



from c_kepler import _kepler as ck
print('imported kepler')
import helper_functions_wrapper as hlpw

import radvel as rv

from scipy.stats import loguniform, beta

## Constants ##

pi    = 3.141592653589793
e     = 2.718281828459045
M_sun = 1.988409870698051e+33
M_jup = 1.8981245973360504e+30


G =  6.674299999999999e-08
au = 14959787070000.0
pc = 3.085677581467192e+18 # 1 parsec in cm

parallax = 18.62 #mas
d_star = (1/parallax)*pc*1000

# @profile
def make_arrays(m_star, a_lim, m_lim, grid_num, num_points):
    
    np.random.seed(0)
    tp = 0
    a_min = a_lim[0]
    a_max = a_lim[1]
    m_min = m_lim[0]
    m_max = m_lim[1]

    
    # These semimajor axes are distances between the planet and the barycenter of the system. The star is on its own orbit, which we will get later.
    a_list = loguniform.rvs(a_min, a_max, size=num_points)
    m_list = loguniform.rvs(m_min, m_max, size=num_points)
    
    # Match up a_list and m_list and get the period for each pair (in days).
    per_list = hlpw.P(a_list, m_star+m_list*(M_jup/M_sun) )
    
    # Eccentricities drawn from a beta distribution. First create a randomized list of probabilities between 0 and 1, then use those probabilities to sample the inverse-CDF of the beta distribution. In essence, for a given cumulative probability, what is the eccentricity? Eccentricities favored by the pdf are more likely to be picked, as desired. I am using (a,b) = (0.867, 3.03) according to Winn & Fabrycky (2014)
    e_list = beta(0.867, 3.03).ppf(np.random.uniform(0, 0.99, num_points))


    cosi_list = np.random.uniform(0, 1, num_points)
    i_list = np.arccos(cosi_list)
    sini_list = np.sqrt(1-cosi_list**2)
    
    
    # Mean anomaly, uniformly distributed. Use this to solve for True anomaly.
    M_anom_list = np.random.uniform(0, 2*pi, num_points)
    # Solving Kepler now to use in RV posterior
    E_anom_list = ck.kepler_array(M_anom_list, e_list)
    T_anom_list = 2*np.arctan(np.sqrt((1+e_list)/(1-e_list)) * np.tan(E_anom_list/2))

    # Arguments of peri, uniformly distributed
    om_list = np.random.uniform(0, 2*pi, num_points)
    
    # Longitudes of ascending node, uniformly distributed
    Om_list = np.random.uniform(0, 2*pi, num_points)

           
    # Breaking up the (a, M) parameter space into grid_num x grid_num
    a_bins = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
    m_bins = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)
    
    a_inds = np.digitize(a_list, bins = a_bins)
    m_inds = np.digitize(m_list, bins = m_bins)
    
    
    
    return a_list, m_list, per_list, e_list, i_list, om_list, M_anom_list, E_anom_list, T_anom_list, a_inds, m_inds

@profile
def rv_post(gammadot, gammadot_err, gammaddot, gammaddot_err, m_star, a_list, m_list, per_list, e_list, i_list, om_list, E_anom_list, num_points, grid_num, a_inds, m_inds):
    
    m_tot_list = (m_star+m_list*(M_jup/M_sun))

    gammadot_list  = np.empty(num_points)
    gammaddot_list = np.empty(num_points)

    # gammadot_list, gammaddot_list = hlpw.gamma(a_list, m_list, per_list, e_list, i_list, om_list, E_anom_list)
    # gammadot_list, gammaddot_list = hlpw.gamma_direct(a_list, m_list, per_list, e_list, i_list, om_list, E_anom_list)
    gammadot_list, gammaddot_list = hlpw.gamma_direct_FAST(a_list, m_list, per_list, e_list, i_list, om_list, E_anom_list)
              
    rv_bounds_memview = hlpw.rv_post_dense_loop(gammadot, gammadot_err, gammaddot, gammaddot_err, gammadot_list, gammaddot_list, a_inds, m_inds, grid_num)


    # Weird situation with a single NaN showing up here. It may have come from the Kepler solver somehow, but I just set all NaNs to 0
    # np.nan_to_num(rv_bounds_array, copy=False)
    rv_bounds_array = np.array(rv_bounds_memview)

    rv_bounds_array = rv_bounds_array/rv_bounds_array.sum()
    
    return rv_bounds_array


def astro_post(delta_mu, delta_mu_err, m_star, a_list, m_list, per_list, e_list, i_list, om_list, T_anom_list, num_points, grid_num, a_inds, m_inds, t_num):
    
    # hip_times  = np.array([2447837.75, 2449065.15])
    # gaia_times = np.array([2456863.5, 2457531.5])
    #
    # hip_baseline = hip_times[1] - hip_times[0]
    # gaia_baseline = gaia_times[1] - gaia_times[0]
    #
    # baseline = ((gaia_times[0] + 0.5*gaia_baseline) - (hip_times[0] + 0.5*hip_baseline))
    #
    #
    # ##########
    # # Experimental: expanding arrays here for cleanliness
    # per_array = per_list[:, None, None]
    # e_array = e_list[:, None, None]
    # T_anom_array = T_anom_list[:, None, None]
    # ##########
    #
    #
    # # Beginning and end points of Hipparcos and Gaia epochs
    # time_endpoints = np.array([[hip_times[0], gaia_times[0]], [hip_times[1], gaia_times[1]]])
    #
    # # Walk through mean anomaly over the HIP/Gaia missions. The mean anomaly array (final shape (per_num, 2, t_num)) contains a mean anomaly for n = (t_num) points along every one of the (per_num) periods, for both Hipparcos and Gaia. per_array and the hip and gaia times are all in days.
    # M_anom_progression = (2*np.pi/per_array)*(np.linspace(time_endpoints[0], time_endpoints[1], t_num) - hip_times[0])[None, :, :]
    # print(M_anom_progression.shape)
    # dfdf
    # M_anom_progression = np.swapaxes(M_anom_progression, 1,2)
    #
    #
    # # Expanding e_array and T_anom_array here to be compatible with the Kepler function
    # e_array = np.repeat(e_array, 2, axis=1)
    # e_array = np.repeat(e_array, 100, axis=2)
    # T_anom_array = np.repeat(T_anom_array, 2, axis=1)
    # T_anom_array = np.repeat(T_anom_array, 100, axis=2)

    x = hlpw.astro_post_dense_loop(delta_mu, delta_mu_err, m_star, a_list, m_list, per_list, e_list, i_list, om_list, T_anom_list, num_points, grid_num, a_inds, m_inds, t_num)
    
    prob_list = []
    astro_bounds_array = np.zeros((grid_num, grid_num))
    
    for i in range(num_points):

        M_2d = M_anom_progression[i]
        e_arr = e_array[i]
        T = T_anom_array[i]
        
        a = a_list[i]
        m = m_list[i]
        per = per_list[i]
        e = e_list[i]
        
        Om = 0 # Om_list[i], remember to verify that this is unnecessary
        om = om_list[i]
        inc = i_list[i]
        
        E_prog_list = []
        for j in range(2):
            
            M_1d = M_2d[j]
            
            # Much faster than hlp.kepler(M_1d, e_arr[0])
            # E_prog = rv._kepler.kepler_array(M_1d, e)
            E_prog = ck.kepler_array(M_1d, e)

            
            E_prog_list.append(E_prog)
        
        # E_prog_list is a (2,100) array with an eccentric anomaly for each of 100 times in both the Hip and Gaia eras.
        E_prog_list = np.array(E_prog_list)

        # Convert the eccentric anomalies into true anomalies and add on T_anom_array, the random starting true anomalies. Transpose so this has shape (t_num, 2). The length-2 axis holds Hipparcos vs. Gaia true anomalies.
        T_prog =  (T + 2*np.arctan( np.sqrt((1+e_arr)/(1-e_arr)) * np.tan(E_prog_list/2))).T
    
    
        # rot_vec = np.array([np.cos(Om)*np.cos(om+T_prog) - np.sin(Om)*np.sin(om+T_prog)*np.cos(inc),
        #                     np.sin(Om)*np.cos(om+T_prog) + np.cos(Om)*np.sin(om+T_prog)*np.cos(inc),
        #                                                              (np.sin(om+T_prog)*np.sin(inc))])
        
        rot_matrix = hlpw.rot_matrix(inc, om, Om)
    
        ######################## Angular Positions ##########################
        # At each of these true anomalies, we need to know the angular separation from barycenter! So first use the r equation on these T_anoms to get r_planet
        r_pl = hlpw.r(T_prog, a*au, e)
        
        # Star physical separation
        r_star = r_pl*((m*M_jup)/(m_star*M_sun))
        
        # r_unit_vec is a stack of arrays that look like [cos(T), sin(T), 0]. It starts with shape (3, 100, 2), but we need it to be 
        # (100, 2, 3, 1) for matrix multiplication, so we move 3 to the end and add another axis. I add negative b/c star is on opposite side of barycenter from planet. Shouldn't actually affect calculation though.
        r_unit_vec = -np.array([np.cos(T_prog), np.sin(T_prog), np.zeros((100,2))])
        r_unit_vec = np.moveaxis(r_unit_vec, 0, 2)[..., None]

        # rot_matrix has shape (3,3), and r_unit_vec is (100, 2, 3, 1), so matrix multuplication ignores (100, 2) and does (3,3)x(3,1) = (3,1). Then we are done with the length-1 axis so we 'squeeze' it out. Finally, move 3 back to the front.
        # rot_vec is the direction of the star relative to the system barycenter on the sky plane. It has norm <= 1.
        rot_vec = np.matmul(rot_matrix, r_unit_vec).squeeze()
        rot_vec = np.moveaxis(rot_vec, 2, 0)

        
        
        # r_star (magnitude) * rot_vec (direction) gives a true vector. Then take only x- and y-components (sky plane). Average over all of the times throughout the epoch to get an average position. Finally, divide by distance to star to convert physical position to angular position.
        # Shape is (2,2): (x&y components, Hip/Gaia)
        ang_pos = ((r_star*rot_vec)[:-1] / d_star * (206265*1e3)).mean(axis= -2) # (cm / cm)*206265*1e3 = mas



        ########################## Angular Velocities ###########################
        
        
        v_vec_pl = hlpw.v_vec(a, per, e, T_prog)
        
        # Stellar velocity is related to planet through their masses. Also they are in opposite directions, so add negative, which shouldn't affect the final answer.

        v_vec_star = -v_vec_pl * ((m*M_jup)/(m_star*M_sun))
        v_vec_star = np.moveaxis(v_vec_star, 0, 2)[..., None]
        
        

        rotated_v_vec = np.matmul(rot_matrix, v_vec_star).squeeze()
        rotated_v_vec = np.moveaxis(rotated_v_vec, 2, 0)
        
        # pm anom values for both Hipparcos and Gaia
        delta_mu_both = ((rotated_v_vec)[:-1] / d_star * (206265*1e3*3.15e7)).mean(axis= -2) # (cm/s / cm)*206265*1e3*3.15e7 = mas/yr
        
        
        # pm_anom_both has shape (2,2): (x&y components, Hip/Gaia)
        # We want only Gaia's x&y components, so we take the whole first axis, and only the second element of the second axis

        delta_mu_gaia = delta_mu_both[:, 1]
        
        
        # Subtract along the second axis, ie, subtract the Gaia and Hipparcos POSITIONS.
        # pm_anom_hg has shape (2), where the 2 components are sin and cos of the hg proper motion anomaly vector.
        # Units are mas/yr
        delta_mu_hg = (ang_pos[:, 1]-ang_pos[:, 0])/(baseline/365)

    
        # Subtract the proper motion anomalies and take the norm along the first axis, which corresponds to the sin and cos components.
        delta_mu_model = np.linalg.norm((delta_mu_gaia-delta_mu_hg), axis = 0)
    
        chi_square = ((delta_mu_model - delta_mu)/delta_mu_err)**2
        
        prob = np.exp(-chi_square/2)
        
        if int(i%(int(num_points/100))) == 0:
            print(int(i / (num_points/100)), '% ')#, prob)
        
        # Placing each probability into the array
        # In case we want 1-sigma bounds for RVs only
        a_i = a_inds[i]
        m_i = m_inds[i]
        astro_bounds_array[m_i, a_i] += prob
        
        
        # prob_list.append(prob) will need this to get the total posterior
    
    # prob_list_astro = np.array(prob_list)
    
    
    # print(np.isnan(self.prob_list_astro).sum())
    
    # np.nan_to_num(self.prob_list_astro, copy=False)

    return astro_bounds_array

    
if __name__ == "__main__":

    gammadot      = 0.114
    gammadot_err  = 0.006
    gammaddot     = -6e-5
    gammaddot_err = 1.9e-5


    m_star = 0.807
    a_lim = (1.9, 5e1)
    m_lim = (1.5, 2e2)
    grid_num = 30
    num_points = int(1e4)
    t_num = 100
    tick_num = 6
    tick_size = 30


    a_list, m_list, per_list, e_list, i_list, om_list, M_anom_list, a_inds, m_inds, other_per_list = make_arrays(m_star, a_lim, m_lim, grid_num, num_points, t_num)


    post_rv = rv_post(gammadot, gammadot_err, gammaddot, gammaddot_err, m_star, a_list, m_list, per_list, e_list, i_list, om_list, M_anom_list, num_points, grid_num, a_inds, m_inds)


    import matplotlib.pyplot as plt
    import matplotlib.patches as ptch


    # # The priors for minimum period and planet mass. min_per is 4xbaseline because we see ~no curvature yet.
    rv_baseline = 430.2527364352718
    max_rv = 40.0581900021484
    min_per = 4*rv_baseline

    # While the above a_list and m_list are the random samples, these are log-uniform lists for plotting.
    a_list = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), grid_num)
    m_list = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), grid_num)

    min_m = rv.utils.Msini(max_rv, min_per, m_star, e=0, Msini_units='jupiter')
    min_a = rv.utils.semi_major_axis(min_per, (m_star + min_m*(M_jup/M_sun)))


    min_index_m = hlpw.value2index(min_m, (0, grid_num-1), m_lim)
    min_index_a = hlpw.value2index(min_a, (0, grid_num-1), a_lim)


    t_contours_rv = hlp.contour_levels(post_rv, [1,2])

    bounds = hlp.bounds_1D(post_rv, [m_lim, a_lim], interp_num = 1e4)
    print('a_lim, m_lim = ', bounds[0], bounds[1])


    fig, ax = plt.subplots(figsize=(12,12))

    post_rv_cont = ax.contourf(post_rv, t_contours_rv, cmap='Greens', extend='max', alpha=0.5)

    mass_rect = ptch.Rectangle((0, 0), grid_num-1, min_index_m, color='gray', alpha=1.0)
    a_rect = ptch.Rectangle((0, 0), min_index_a, grid_num-1, color='gray', alpha=1.0)

    ax.add_patch(mass_rect)
    ax.add_patch(a_rect)

    tick_array = np.linspace(0, grid_num-1, tick_num).astype(int)

    plt.xticks(tick_array, [np.round(a_list[i], 1) for i in tick_array], size=tick_size)
    plt.yticks(tick_array, [np.round(m_list[i], 1) for i in tick_array ], size=tick_size)


    fig.tight_layout()
    # fig.savefig('5thCompConstraints_RV_astr.png')
    plt.show()
    
    
    
    
    
    
    
    