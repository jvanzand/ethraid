import os
from astropy.time import Time

__version__='2.5.2' # Remember to remove random seed

_ROOT = os.path.abspath(os.path.dirname(__file__))

# CONSTANTS
# Mass of the Sun in grams
M_sun = 1.988409870698051e+33
# Mass of Jupiter in grams
M_jup = 1.8981245973360504e+30
# 1 parsec in AU
pc_in_au = 206264.80624548031 # (c.pc.cgs/c.au.cgs).value
# Conversion factor between solar masses and Jupiter masses
Ms2Mj = M_sun/M_jup

# Astrometry times
#https://www.cosmos.esa.int/web/hipparcos/catalogue-summary
hip_times  = [Time(1989.85, format='decimalyear').jd, Time(1993.21, format='decimalyear').jd]
hip_epoch = 1991.25

#https://www.cosmos.esa.int/web/gaia/earlydr3
gaia_times = [Time('2014-07-25', format='isot').jd, Time('2017-05-28', format='isot').jd]
gaia_epoch = 2016.0