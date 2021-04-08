import matplotlib.pyplot as plt
import astropy.constants as c
import numpy as np
import helper_functions as hlp

import radvel as rv
from c_kepler import _kepler as ck

from scipy.stats import loguniform, beta

from constants import *


"""
This module is an experimental version of giant_class.py. 
I am trying completely randomized points rather than moving through arrays.
Edit: this version is now OFFICIAL
"""
class Giant(object):
    
    def __init__(self, star_name, m_star, gamma_dot, gamma_dot_err, 
                gamma_dotdot, gamma_dotdot_err, parallax, pm_anom_data, pm_anom_data_err):
        
        self.star_name = star_name
        self.m_star = m_star
        self.gamma_dot = gamma_dot
        self.gamma_dot_err = gamma_dot_err
        self.gamma_dotdot = gamma_dotdot
        self.gamma_dotdot_err = gamma_dotdot_err
        self.parallax = parallax # mas
        self.d_star = (1/parallax)*c.pc.cgs.value*1000 # cm
        self.pm_anom_data = pm_anom_data
        self.pm_anom_data_err = pm_anom_data_err
        
        
    
    def joint_post_rv(self, a_lim = (1.9, 1e2), m_lim = (1.5, 1e2), grid_num = 100, t_num = 100, e = 0, tp = 0, om = 0):

        a_min, a_max = a_lim
        m_min, m_max = m_lim
    
        a_list = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
        m_list = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)

        prob_array = np.zeros((len(m_list), len(a_list)))
    
    
        for j in range(len(a_list)):
            for i in range(len(m_list)): # I'm assuming 90 inclination for simplicity, so these will be our "msini" values. They are still in Jupiter masses
            
                a = a_list[j]
                m = m_list[i]
            
                per = hlp.P(a, (self.m_star+m*(c.M_jup/c.M_sun).value))
                k = rv.utils.semi_amplitude(m, per, (self.m_star+m*(c.M_jup/c.M_sun).value), e, Msini_units='jupiter')

                orbel = [per, tp, e, om, k]


                t_list = np.linspace(0, per, t_num) # List of times spanning one full orbit. Dates are arbitrary so I start at jd=0 for ease

                rv_list = rv.kepler.rv_drive(t_list, orbel)

                gamma_dot_array = hlp.deriv(rv_list, t_list)[1:-1]
                gamma_dotdot_array = hlp.deriv(hlp.deriv(rv_list, t_list), t_list[1:-1])
        
                dot_term = ((gamma_dot_array-self.gamma_dot)/(self.gamma_dot_err))**2
                dotdot_term = ((gamma_dotdot_array-self.gamma_dotdot)/(self.gamma_dotdot_err))**2
                probs = np.e**(-(dot_term + dotdot_term)/2)
                
                prob = np.sum(probs) # By summing across the orbit, we effectively marginalize over time of conjunction.

                prob_array[i,j] = prob

        prob_array /= prob_array.sum()
    
        return prob_array
        
    def joint_post_astro(self, a_lim = (1.9, 1e2), m_lim = (1.5, 1e2), grid_num = 102):
        
        a_min, a_max = a_lim
        m_min, m_max = m_lim
    
        a_list = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
        m_list = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)

        prob_array = np.zeros((len(m_list), len(a_list)))

        prob_sum = 0
        for j in range(len(a_list)):
            for i in range(len(m_list)):
            
                a = a_list[j]
                m = m_list[i]
                
            
                # The planet and the star have the same period
                per = hlp.P(a, (self.m_star+m*(c.M_jup/c.M_sun).value))
                per_sec = per*(24*3600)
            
                # Astrometry gives us the STELLAR proper motions, so we need to simulate those. We start by getting the separation of the star from the system COM
                a_star_cgs = (a*c.au.cgs.value)*((m*c.M_jup.cgs.value)/(self.m_star*c.M_sun.cgs.value))
                v_cgs = (2*np.pi*a_star_cgs)/per_sec

            
                # We need to convert a_star_cgs and v_cgs to angular values (mas and mas/yr)
                ang_sep = (a_star_cgs/self.d_star)*(206265*1e3)
                pm_anom = (v_cgs/self.d_star)*(206265*1e3)*(3.15e7)
                
            
                # Find the true anomaly of the planet at each epoch. Since Tim Brandt's catalog is dealing with AVERAGE positions and velocities over the Hip/Gaia missions, that's what I calculate here. I chose t = 0 at the beginning of the Hipparcos mission. 
            
                nu_hip = (2*np.pi/per)*(np.linspace(hip_times[0], hip_times[1], 1000) - hip_times[0])
                nu_gaia = (2*np.pi/per)*(np.linspace(gaia_times[0], gaia_times[1], 1000) - hip_times[0])
            
            
                # Positions in mas of the star at gaia/hip times to calculate HG velocity
                # I'm treating the barycenter as (0,0) to get the difference in position
                # Calculate many positions, then average over the x- and y- components to get 1 position

                pos_hip = ang_sep*np.array((np.cos(nu_hip), np.sin(nu_hip))).mean(1)
                pos_gaia = ang_sep*np.array((np.cos(nu_gaia), np.sin(nu_gaia))).mean(1)
            
                
                # if i == 5 and j == 0:
                #     print(m*c.M_jup.cgs.value)
                #     print(a_star_cgs)
                   
                # For a circular orbit, the pm_anomaly components depend directly on nu
            
                # v_hip = pm_anom*np.array((-np.sin(nu_hip), np.cos(nu_hip))).mean(1)
                v_gaia = pm_anom*np.array((-np.sin(nu_gaia), np.cos(nu_gaia))).mean(1)
            
                v_hg = (pos_gaia-pos_hip)/(baseline/365)
            
                # Now we have the acceleration of the host star in mas/yr/day
                # Switching to delta_pm, we remove the division by baseline
                pm_anom_model = np.linalg.norm((v_gaia-v_hg))

                chi_square = ((pm_anom_model - self.pm_anom_data)/self.pm_anom_data_err)**2
            
                prob = np.exp(-chi_square/2)
        
                prob_array[i,j] = prob
                
    
        prob_array /= prob_array.sum()
        return prob_array
    
    def make_arrays(self, a_lim = (1.9, 1e2), m_lim = (1.5, 1e2), grid_num = 100, num_points = int(3e5), t_num = 100, plot_num = 30, e = 0, tp = 0, om = 0):
        
        

        self.tp = tp
        np.random.seed(0)
        
        a_min, a_max = a_lim
        m_min, m_max = m_lim
        
        self.num_points = num_points
        self.grid_num = grid_num
        self.plot_num = plot_num
        # These semimajor axes are distances between the planet and the barycenter of the system. The star is on its own orbit, which we will get later.
        self.a_list = loguniform.rvs(a_min, a_max, size=num_points)
        self.m_list= loguniform.rvs(m_min, m_max, size=num_points)
        
        # Match up a_list and m_list and get the period for each pair.
        self.per_list = hlp.P(self.a_list, (self.m_star+self.m_list*(c.M_jup/c.M_sun).value))
        
        # Eccentricities drawn from a beta distribution. First create a randomized list of probabilities between 0 and 1, then use those probabilities to sample the inverse-CDF of the beta distribution. I am using (a,b) = (0.867, 3.03) according to Winn & Fabrycky (2014)
        self.e_list = beta(0.867, 3.03).ppf(np.random.uniform(0, 0.99, num_points))

    
        self.cosi_list = np.random.uniform(0, 1, num_points)
        self.i_list = np.arccos(self.cosi_list)
        self.sini_list = np.sqrt(1-self.cosi_list**2)
        
        
        # Mean anomaly, uniformly distributed. Use this to solve for True anomaly.
        self.M_anom_list = np.random.uniform(0, 2*np.pi, num_points)
        # Solving Kepler now to use in RV posterior
        self.E_anom_list = ck.kepler_array(self.M_anom_list, self.e_list)
        self.T_anom_list = 2*np.arctan(np.sqrt((1+self.e_list)/(1-self.e_list)) * np.tan(self.E_anom_list/2))

        # Arguments of peri, uniformly distributed
        self.om_list = np.random.uniform(0, 2*np.pi, num_points)
        
        # Longitudes of ascending node, uniformly distributed
        self.Om_list = np.random.uniform(0, 2*np.pi, num_points)

               
        # Breaking up the (a, M) parameter space into grid_num x grid_num
        a_bins = np.logspace(np.log10(a_min), np.log10(a_max), grid_num)
        m_bins = np.logspace(np.log10(m_min), np.log10(m_max), grid_num)
        
        self.a_inds = np.digitize(self.a_list, bins = a_bins)
        self.m_inds = np.digitize(self.m_list, bins = m_bins)
        
        # Breaking up the (a, M) parameter space into plot_num x plot_num
        a_bins_plot = np.logspace(np.log10(a_min), np.log10(a_max), plot_num)
        m_bins_plot = np.logspace(np.log10(m_min), np.log10(m_max), plot_num)
        
        self.a_inds_plot = np.digitize(self.a_list, bins = a_bins_plot)
        self.m_inds_plot = np.digitize(self.m_list, bins = m_bins_plot)
        
    
        # self.msini_list = self.m_list*self.sini_list
        
        
        return
    
    def rv_post(self):
        
        m_tot_list = (self.m_star+self.m_list*(c.M_jup/c.M_sun).value)
        
        # Get the first and second derivatives of the RV curve
        gamma_dot, gamma_dotdot = hlp.gamma(self.a_list, self.m_list, self.per_list, self.e_list, self.i_list, self.om_list, self.M_anom_list)
        
        dot_term_list = ((gamma_dot-self.gamma_dot)/(self.gamma_dot_err))**2
        dotdot_term_list = ((gamma_dotdot-self.gamma_dotdot)/(self.gamma_dotdot_err))**2
        
        # Why three arrays to store probabilities? self.prob_list_rv is a 1D list that will be multiplied with its astrometry counterpart to get a the total (rv+astro) marginalized posterior. rv_bounds_array is a 2D array of the rv probabilities, binned according to a and M. It will be used to find the 1-sigma bounds on a an M. rv_plot_array is similar, but it will be used for plotting the 2D probability surface instead of finding bounds. The only difference between the last two is the number of points used to create the 2D grid.
        
        self.prob_list_rv = np.exp(-(dot_term_list + dotdot_term_list)/2)
        # Weird situation with a single NaN showing up here. It may have come from the Kepler solver somehow, but I just set all NaNs to 0
        np.nan_to_num(self.prob_list_rv, copy=False)
        
        
        rv_bounds_array = np.zeros((self.grid_num, self.grid_num))
        rv_plot_array   = np.zeros((self.plot_num, self.plot_num))
        
        for i in range(self.num_points):
          # In case we want 1-sigma bounds for RVs only
          a_i = self.a_inds[i]
          m_i = self.m_inds[i]
          rv_bounds_array[m_i, a_i] += self.prob_list_rv[i]
          
          # For plotting
          a_i_plot = self.a_inds_plot[i]
          m_i_plot = self.m_inds_plot[i]
          rv_plot_array[m_i_plot, a_i_plot] += self.prob_list_rv[i]
          
          if int(i%(int(self.num_points/100))) == 0:
              print(int(i / (self.num_points/100)), '% ', self.prob_list_rv[i])

        self.rv_bounds_array = rv_bounds_array/rv_bounds_array.sum()
        self.rv_plot_array   = rv_plot_array/rv_plot_array.sum()
        return

    def astro_post_new(self):

        ##########
        # Experimental: expanding arrays here for cleanliness
        per_array = self.per_list[:, None, None]
        e_array = self.e_list[:, None, None]
        T_anom_array = self.T_anom_list[:, None, None]

        ##########
        
        # Beginning and end points of Hipparcos and Gaia epochs
        time_endpoints = np.array([[hip_times[0], gaia_times[0]], [hip_times[1], gaia_times[1]]])

        # Walk through mean anomaly over the HIP/Gaia missions. The mean anomaly array (shape (per_num, t_num, 2)) contains a mean anomaly for n = (t_num) points along every one of the (per_num) periods, for both Hipparcos and Gaia.
        

        M_anom_progression = (2*np.pi/self.per_list)*(np.linspace(time_endpoints[0], time_endpoints[1], 100) - hip_times[0])[..., None]
        

        # Have to jump through some hoops to calculate E anomaly because the C Kepler solver can only take 1-D arrays. See if I can ask Erik about changing this. Maybe not, because he said C trades versatility for speed, ie, you can only supply a function with a single type of input. And if that's true, then a 1-D array is better than a 3-D or whatever I would use instead.
        E_anom_progression = []
        for i in range(2):
            epoch_progression = []
            for j in range(self.num_points):
                # M is a list of mean anomalies for either Hip or Gaia and for one set of random parameters
                M_sublist = M_anom_progression[:, i, j]
                
                # e is a list of one single value, the random eccentricity corresponding to the random parameters of M
                e_sublist = np.repeat(self.e_list[i], 100)
                
                
                E_sublist = ck.kepler_array(M_sublist, e_sublist)

                epoch_progression.append(E_sublist)
                
                if int(j%(int(self.num_points/100))) == 0:
                    print(int(j / (self.num_points/100)), '% ')
            
            E_anom_progression.append(epoch_progression)
        

        E_anom_progression = np.swapaxes(np.array(E_anom_progression), 1,2)
        

        
        
        # Convert the eccentric anomalies into true anomalies and add on T_anom_array, the random starting true anomalies. Transpose so this has shape (2, t_num, per_num) and can multiply other arrays with shape (per_num). The length-2 axis holds Hipparcos vs. Gaia true anomalies.
        # T_prog =  (T + 2*np.arctan( np.sqrt((1+e_arr)/(1-e_arr)) * np.tan(E_prog/2))).T
        
        T_prog = (self.T_anom_list + 2*np.arctan( np.sqrt((1+self.e_list)/(1-self.e_list)) * np.tan(E_anom_progression/2)))
    
        
        # rot_vec has shape (3, 2, 100, num_points). 3 spatial components, 2 epochs, 100 t points, and num_points of random angle combinations.
        rot_vec = np.array([np.cos(self.Om_list)*np.cos(self.om_list+T_prog) - np.sin(self.Om_list)*np.sin(self.om_list+T_prog)*np.cos(self.i_list),
                            np.sin(self.Om_list)*np.cos(self.om_list+T_prog) + np.cos(self.Om_list)*np.sin(self.om_list+T_prog)*np.cos(self.i_list),
                                                                                                   (np.sin(self.om_list+T_prog)*np.sin(self.i_list))])

    
        # At each of these true anomalies, we need to know the angular separation! So use the r equation on these nu's to get r_planet
    
        r_pl = hlp.r(T_prog, self.a_list*c.au.cgs.value, self.e_list)

        r_star = r_pl*((self.m_list*c.M_jup.cgs.value)/(self.m_star*c.M_sun.cgs.value))
        
        
        # (r_star*rot_vec) has shape (3, 2, 100, num_points). We only want the x,y components of the vector, so omit the z component. Then take the norm over the x and y.
        
        # proj_dist = np.linalg.norm((r_star*rot_vec)[:-1], axis = 0)
        # ang_sep = (proj_dist/self.d_star)*(206265*1e3)

        
        # Angular positions varying along t_num = 100. I average over the epoch times to get an average position on the sky.
        ang_pos = ((r_star*rot_vec)[:-1] / self.d_star * (206265*1e3)).mean(axis= -2)

        
        # These lines are analogous to the 4 above, except now instead of the magnitude of angular position (ang_sep), we want the magnitude of angular velocity (pm_mag).
        r_dot_pl = hlp.r_dot(T_prog, self.a_list*c.au.cgs.value, self.per_list*(24*3600), self.e_list)
    
        r_dot_star = r_dot_pl * ((self.m_list*c.M_jup.cgs.value)/(self.m_star*c.M_sun.cgs.value))
    
        proj_v = np.linalg.norm((r_dot_star*rot_vec)[:-1], axis = 0)
    
        pm_mag = (proj_v/self.d_star)*(206265*1e3)*(3.15e7)

    

        # Average over t_num, the different times in the Hip and Gaia eras. This leaves behind a (2,2, num_points) array, which are sin and cos components for Hip and Gaia for each set of random parameters.
        
        
        
        # ang_pos = (ang_sep*np.array((np.cos(T_prog), np.sin(T_prog)))).mean(-2)
        

        # Subtract along the second axis, ie, subtract the Gaia and Hipparcos positions.
        # pm_anom_hg has shape (2, num_points), where the 2 components are sin and cos of the hg proper motion anomaly vector.
        pm_anom_hg = (ang_pos[:, 1]-ang_pos[:, 0])/(baseline/365)

    
        # We want only the gaia mean proper motion anomaly, so we take only the Gaia true anomalies. Then we multiply these by only the Gaia proper motion magnitudes.

        
        # Take only the Gaia pm_mag values, and multiply these by the sin and cos components of the Gaia true anomalies.
        # This gives a shape (2, num_points) array, whose components are the sin and cos of the pm_anom in the Gaia epoch.

        pm_anom_gaia = (pm_mag[1]*np.array((-np.sin(T_prog[1]), np.cos(T_prog[1])))).mean(-2)

    
        # Subtract the proper motion anomalies and take the norm along the first axis, which corresponds to the sin and cos components.
        pm_anom_model = np.linalg.norm((pm_anom_gaia-pm_anom_hg), axis = 0)
    
        chi_square = ((pm_anom_model - self.pm_anom_data)/self.pm_anom_data_err)**2
        
        self.prob_list_astro = np.exp(-chi_square/2)
        
        
        prob_array = np.zeros((self.grid_num, self.grid_num))
        plot_array = np.zeros((self.plot_num, self.plot_num))
        for i in range(self.num_points):
            
          a_i = self.a_inds[i]
          m_i = self.m_inds[i]
          
          prob_array[m_i, a_i] += self.prob_list_astro[i]

        return prob_array/prob_array.sum()
        
        
    
    def astro_post(self):

        ##########
        # Experimental: expanding arrays here for cleanliness
        per_array = self.per_list[:, None, None]
        e_array = self.e_list[:, None, None]
        T_anom_array = self.T_anom_list[:, None, None]
        ##########
        
        
        # Beginning and end points of Hipparcos and Gaia epochs
        time_endpoints = np.array([[hip_times[0], gaia_times[0]], [hip_times[1], gaia_times[1]]])

        # Walk through mean anomaly over the HIP/Gaia missions. The mean anomaly array (final shape (per_num, 2, t_num)) contains a mean anomaly for n = (t_num) points along every one of the (per_num) periods, for both Hipparcos and Gaia.
        M_anom_progression = (2*np.pi/per_array)*(np.linspace(time_endpoints[0], time_endpoints[1], 100) - hip_times[0])[None, :, :]
        M_anom_progression = np.swapaxes(M_anom_progression, 1,2)
        
        
        # Expanding e_array and T_anom_array here to be compatible with the Kepler function
        e_array = np.repeat(e_array, 2, axis=1)
        e_array = np.repeat(e_array, 100, axis=2)
        T_anom_array = np.repeat(T_anom_array, 2, axis=1)
        T_anom_array = np.repeat(T_anom_array, 100, axis=2)

        
        prob_list = []
        astro_bounds_array = np.zeros((self.grid_num, self.grid_num))
        astro_plot_array   = np.zeros((self.plot_num, self.plot_num))
        for i in range(self.num_points):

            M_2d = M_anom_progression[i]
            e_arr = e_array[i]
            T = T_anom_array[i]
            
            a = self.a_list[i]
            m = self.m_list[i]
            per = self.per_list[i]
            e = self.e_list[i]
            
            Om = self.Om_list[i]
            om = self.om_list[i]
            inc = self.i_list[i]
            

            E_prog_list = []
            for j in range(2):
                
                M_1d = M_2d[j]
                E_prog = rv._kepler.kepler_array(M_1d, e)
                # E_prog = hlp.kepler(M_1d, e_arr[0])

                
                E_prog_list.append(E_prog)
            
            # E_prog_list is a (2,100) array with an eccentric anomaly for each of 100 times in both the Hip and Gaia eras.
            E_prog_list = np.array(E_prog_list)
            

            # Convert the eccentric anomalies into true anomalies and add on T_anom_array, the random starting true anomalies. Transpose so this has shape (t_num, 2). The length-2 axis holds Hipparcos vs. Gaia true anomalies.
            T_prog =  (T + 2*np.arctan( np.sqrt((1+e_arr)/(1-e_arr)) * np.tan(E_prog_list/2))).T
        
        
            rot_vec = np.array([np.cos(Om)*np.cos(om+T_prog) - np.sin(Om)*np.sin(om+T_prog)*np.cos(inc),
                                np.sin(Om)*np.cos(om+T_prog) + np.cos(Om)*np.sin(om+T_prog)*np.cos(inc),
                                                                         (np.sin(om+T_prog)*np.sin(inc))])

        
            ######################## Angular Positions ##########################
            # At each of these true anomalies, we need to know the angular separation from barycenter! So first use the r equation on these T_anoms to get r_planet
            r_pl = hlp.r(T_prog, a*c.au.cgs.value, e)

            # Star physical separation
            r_star = r_pl*((m*c.M_jup.cgs.value)/(self.m_star*c.M_sun.cgs.value))
    
            # r_star * rot_vec gives 1) distance with 2) direction. Then take only x- and y-components (sky plane). Average over all of the times throughout the epoch to get an average position. Finally, divide by distance to star to convert physical position to angular position.
            # Shape is (2,2): (x&y components, Hip/Gaia)

            ang_pos = ((r_star*rot_vec)[:-1] / self.d_star * (206265*1e3)).mean(axis= -2) # (cm / cm)*206265*1e3 = mas

        
            ########################## Angular Velocities ###########################
            
            r_dot_pl = hlp.r_dot(T_prog, a*c.au.cgs.value, per*(24*3600), e)
        
            r_dot_star = r_dot_pl * ((m*c.M_jup.cgs.value)/(self.m_star*c.M_sun.cgs.value))
            
            # pm_anom_both has shape (2,2): (x&y components, Hip/Gaia)
            pm_anom_both = ((r_dot_star*rot_vec)[:-1] / self.d_star * (206265*1e3*3.15e7)).mean(axis= -2) # (cm/s / cm)*206265*1e3*3.15e7 = mas/yr
            
            
            
            # We want only Gaia's x&y components, so we take the second element of the second axis
            pm_anom_gaia = pm_anom_both[:, 1]
            
            
            # Subtract along the second axis, ie, subtract the Gaia and Hipparcos POSITIONS.
            # pm_anom_hg has shape (2), where the 2 components are sin and cos of the hg proper motion anomaly vector.
            # Units are mas/yr
            pm_anom_hg = (ang_pos[:, 1]-ang_pos[:, 0])/(baseline/365)

        
            # Subtract the proper motion anomalies and take the norm along the first axis, which corresponds to the sin and cos components.
            pm_anom_model = np.linalg.norm((pm_anom_gaia-pm_anom_hg), axis = 0)
        
            chi_square = ((pm_anom_model - self.pm_anom_data)/self.pm_anom_data_err)**2
            
            prob = np.exp(-chi_square/2)
            
            if int(i%(int(self.num_points/100))) == 0:
                print(int(i / (self.num_points/100)), '% ', prob)
            
            # Placing each probability into the array
            # In case we want 1-sigma bounds for RVs only
            a_i = self.a_inds[i]
            m_i = self.m_inds[i]
            astro_bounds_array[m_i, a_i] += prob
            
            
            # For plotting
            a_i_plot = self.a_inds_plot[i]
            m_i_plot = self.m_inds_plot[i]
            astro_plot_array[m_i_plot, a_i_plot] += prob
            
            prob_list.append(prob)
    
        self.prob_list_astro = np.array(prob_list)
        
        self.astro_bounds_array = astro_bounds_array
        self.astro_plot_array   = astro_plot_array
        
        print(np.isnan(self.prob_list_astro).sum())
        
        np.nan_to_num(self.prob_list_astro, copy=False)

        return
        

    
    def rv_astro_post(self):
        """
        To get the combined RV/astrometry posterior, we need to multiply the RV/astrometry probabilities at EACH point,
        and THEN marginalize.
        """

        multiplied_probs = np.array(self.prob_list_rv) * np.array(self.prob_list_astro)

        prob_array            = np.zeros((self.grid_num, self.grid_num))
        self.tot_plot_array = np.zeros((self.plot_num, self.plot_num))
        for i in range(self.num_points):
            
            prob = multiplied_probs[i]
            
            a_i = self.a_inds[i]
            m_i = self.m_inds[i]
            prob_array[m_i, a_i] += prob
            
            # For plotting
            a_i_plot = self.a_inds_plot[i]
            m_i_plot = self.m_inds_plot[i]
            self.tot_plot_array[m_i_plot, a_i_plot] += prob
    
        return prob_array/prob_array.sum()
    
if __name__ == "__main__":
    
    grid_num = 50
    plot_num = 30
    num_points = int(1e4)
    tick_num = 6
    tick_size = 30
    label_size = 36
    title_size = 38
    
    a_lim = (1.9, 5e1)
    m_lim = (1.5, 2e2)
    
    
    hd191939 = ['HD191939', 0.807, 0.114, 0.006, -6e-5, 1.9e-5, 18.62, 0.1187, 0.0961]
    
    my_planet = Giant(*hd191939)
    
    my_planet.make_arrays(grid_num = grid_num, num_points = num_points, plot_num = plot_num)
    #
    # # old_astro = my_planet.joint_post_astro()
    my_planet.astro_post()
    post_astro = my_planet.astro_plot_array
    # post_rv = my_planet.rv_post()
    # post_tot = my_planet.rv_astro_post()


    # fig, axs = plt.subplots(1,3)
    fig, ax = plt.subplots(figsize=(15,12))
    ax.imshow(post_astro, origin='lower', cmap='jet')
    
    tick_array = np.linspace(0, plot_num-1, tick_num+1).astype(int)
    
    a_list = np.logspace(np.log10(a_lim[0]), np.log10(a_lim[1]), plot_num)
    m_list = np.logspace(np.log10(m_lim[0]), np.log10(m_lim[1]), plot_num)
    
    plt.xticks(tick_array, [np.round(a_list[i], 1) for i in tick_array], size=tick_size)
    plt.yticks(tick_array, [np.round(m_list[i], 1) for i in tick_array ], size=tick_size)
    
    ax.set_xlabel('Semi-major Axis (au)', size=label_size)
    ax.set_ylabel(r'M ($M_{Jup}$)', size=label_size)
    ax.set_title('Astro Constraints', size=title_size)

    # axs[1].imshow((post_astro + post_rv) / (post_astro + post_rv).sum(), origin='lower', cmap = 'jet')
    # axs[1].imshow(post_rv, origin='lower', cmap = 'jet')
    # axs[2].imshow(post_tot, origin='lower', cmap = 'jet')
    #
    plt.show()



























