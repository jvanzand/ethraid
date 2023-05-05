import os

__version__='2.1.0'

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