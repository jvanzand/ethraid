from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np
import re

import line_profiler
profile = line_profiler.LineProfiler()

####################
# # Extra stuff to profile code using kernprof -l -v file.py after it is compiled
# # Taken from https://stackoverflow.com/questions/28301931/how-to-profile-cython-functions-line-by-line
#
# from Cython.Compiler.Options import get_directive_defaults
# directive_defaults = get_directive_defaults()
# import Cython
# directive_defaults = Cython.Compiler.Options.get_directive_defaults()
#
# directive_defaults['linetrace'] = True
# directive_defaults['binding'] = True

###########################################

# First path is where to put compiled files (.so).
# Second is where to find .pyx files
extensions = [
    Extension("c_kepler._kepler",
             ['c_kepler/_kepler.pyx'],
              include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
    Extension("trends.compiled.helper_functions_general",
             ['trends/helper_functions_general.pyx'],
              include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
    Extension("trends.compiled.helper_functions_rv",
             ['trends/helper_functions_rv.pyx'],
              include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
    Extension("trends.compiled.helper_functions_astro",
             ['trends/helper_functions_astro.pyx'],
              include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')])
]

##############################

# # Specify a build directory to store the .c and .html files separately. Improves readability.
# setup(
#     ext_modules = cythonize(extensions, build_dir = 'trends/compiled', language_level='3', annotate=True)
# )

setup(
    name = "trends",
    packages = ["trends"],
    entry_points = {
        "console_scripts": ['trends = trends.cli:main']
        },
    version = "0.1.0",
    description = "Command line application for trend code",
    author = "Judah Van Zandt",
    author_email = "judahvz@astro.ucla.edu",
    url = "",
    ext_modules = cythonize(extensions, build_dir = 'trends/compiled', language_level='3', annotate=True)
    )
    
    
    
    