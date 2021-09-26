from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np

# import line_profiler
# profile = line_profiler.LineProfiler()

####################
# Extra stuff to profile code using kernprof -l -v file.py after it is compiled
# Taken from https://stackoverflow.com/questions/28301931/how-to-profile-cython-functions-line-by-line

# # from Cython.Compiler.Options import get_directive_defaults
# # directive_defaults = get_directive_defaults()
# import Cython
# directive_defaults = Cython.Compiler.Options.get_directive_defaults()
#
# directive_defaults['linetrace'] = True
# directive_defaults['binding'] = True

###########################################

# First path is where to put build files (.c and .html).
# Second is where to 1) find .pyx files and 2) store .so files.
extensions = [
    Extension("c_kepler._kepler",
             ['c_kepler/_kepler.pyx'],
             include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
    Extension("trends.helper_functions_general",
             ['trends/helper_functions_general.pyx'],
             include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
    Extension("trends.helper_functions_rv",
             ['trends/helper_functions_rv.pyx'],
             include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
    Extension("trends.helper_functions_astro",
             ['trends/helper_functions_astro.pyx'],
             include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')])
    
]

##############################

# Specify a build directory to store the .c and .html files separately. The .so files are what I need to call, so those get put in the same directory as their .pyx file. Improves readability.
setup(
    ext_modules = cythonize(extensions, build_dir = 'trends/compiled', language_level='3', annotate=True)
)