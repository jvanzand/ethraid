import numpy as np
import matplotlib.pyplot as plt
import radvel as rv
import helper_functions as hlp

n = int(1e4)
a = 1 # au
Mp = 1 # M_Jup
per = 365.25 # days
m_star = 1047.5 # M_Jup, ~ 1 M_Sun
i = np.pi/2
e = 0.8
# M_anom_list = np.linspace(0, 2*np.pi, n)
M_anom_list = np.linspace(-np.pi, np.pi, n)
# M_anom_list = M_anom_list*np.array([1,2,3])
# print(M_anom_list)
e_list = np.array([e for i in range(n)])

# E_anom_list = hlp.kepler(M_anom_list, e)
# T_anom_list = 2*np.arctan(np.sqrt((1+e)/(1-e)) * np.tan(E_anom_list/2))

om = 0


### Numerical ###

k = rv.utils.semi_amplitude(Mp*np.sin(i), per, (m_star + Mp)/1047.5, e, Msini_units='jupiter')

tp = 0
orbel = [per, tp, e, om, k]

# t_list = np.linspace(0, per, 1000)
t_list = np.linspace(-per/2, per/2, n)
rv_array = rv.kepler.rv_drive(t_list, orbel)




fig, axs = plt.subplots(2,2)
fig.tight_layout(pad=4)



gd, gdd = hlp.gamma(a, Mp, per, e_list, i, om, M_anom_list)


print(gd.shape)


axs[0,0].plot(M_anom_list, gd)
axs[1,0].plot(t_list[2:-2], hlp.deriv(rv_array, t_list)[1:-1])

axs[0,1].plot(M_anom_list, gdd)
axs[1,1].plot(t_list[2:-2], hlp.deriv(hlp.deriv(rv_array, t_list), t_list[1:-1]))


axs[0,0].set_xlabel(r'$M_{anom}$')
axs[0,0].set_ylabel(r'$\dot{\gamma}$ (m/s/d)')
axs[0,0].set_title(r'$\dot{\gamma}(M)$ (analytic)')

axs[1,0].set_xlabel('time since tp (days)')
axs[1,0].set_ylabel(r'$\dot{\gamma}$ (m/s/d)')
axs[1,0].set_title(r'$\dot{\gamma}(t)$ (numerical)')

axs[0,1].set_xlabel(r'$M_{anom}$')
axs[0,1].set_ylabel(r'$\ddot{\gamma} (m/s/d^{2})$')
axs[0,1].set_title(r'$\ddot{\gamma}(M)$ (analytic)')

axs[1,1].set_xlabel('time since tp (days)')
axs[1,1].set_ylabel(r'$\ddot{\gamma} (m/s/d^{2})$')
axs[1,1].set_title(r'$\ddot{\gamma}(t)$ (numerical)')


plt.show()













