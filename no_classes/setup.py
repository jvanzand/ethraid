from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension
import numpy as np

# import line_profiler
# profile = line_profiler.LineProfiler()

####################

# Extra stuff to profile code using kernprof -l -v file.py after it is compiled
# Taken from https://stackoverflow.com/questions/28301931/how-to-profile-cython-functions-line-by-line

from Cython.Compiler.Options import get_directive_defaults
directive_defaults = get_directive_defaults()

directive_defaults['linetrace'] = True
directive_defaults['binding'] = True

#####################

# Extension(name_of_so_file, [name_of_pyx_file.pyx])

##############################
extensions = [
    Extension("helper_functions_wrapper", ['helper_functions_wrapper.pyx'], include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
    Extension("helper_functions", ['helper_functions.pyx'], include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')])
]
##############################

setup(
    ext_modules = cythonize(extensions, language_level='3', annotate=False)
)