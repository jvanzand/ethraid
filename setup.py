from setuptools import setup, find_packages, Extension
import subprocess
from setuptools.command.build_ext import build_ext as _build_ext
import re

import os
import sys

####################
# ## Extra stuff to profile code using kernprof -l -v file.py after it is compiled
# ## Taken from https://stackoverflow.com/questions/28301931/how-to-profile-cython-functions-line-by-line
#
# import line_profiler
# profile = line_profiler.LineProfiler()
# import numpy as np
# from Cython.Compiler.Options import get_directive_defaults
# directive_defaults = get_directive_defaults()
#
# directive_defaults['linetrace'] = True
# directive_defaults['binding'] = True
#
# extensions = [
#     Extension("ethraid.compiled._kepler",
#              ['ethraid/c_kepler/_kepler.pyx'], include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
#     Extension("ethraid.compiled.helper_functions_general",
#              ['ethraid/helper_functions_general.pyx'], include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
#     Extension("ethraid.compiled.helper_functions_rv",
#              ['ethraid/helper_functions_rv.pyx'], include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
#     Extension("ethraid.compiled.helper_functions_astro",
#              ['ethraid/helper_functions_astro.pyx'], include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
#     Extension("ethraid.compiled.helper_functions_imaging",
#              ['ethraid/helper_functions_imaging.pyx'], include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')])
# ]

###########################################

# First path is where to put compiled files (.so).
# Second is where to find .pyx files
## If profiling, add: include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]
## to each extension

extensions = [
    Extension("ethraid.compiled._kepler",
             ['ethraid/c_kepler/_kepler.pyx']),
    Extension("ethraid.compiled.helper_functions_general",
             ['ethraid/helper_functions_general.pyx']),
    Extension("ethraid.compiled.helper_functions_rv",
             ['ethraid/helper_functions_rv.pyx']),
    Extension("ethraid.compiled.helper_functions_astro",
             ['ethraid/helper_functions_astro.pyx']),
    Extension("ethraid.compiled.helper_functions_imaging",
             ['ethraid/helper_functions_imaging.pyx'])       
]

##############################

try:
    from Cython.Build import cythonize
    import numpy
except ModuleNotFoundError:
    subprocess.run(["pip", "install", "cython==0.29.33"])
    subprocess.run(["pip", "install", "numpy==1.21.6"])
    from Cython.Build import cythonize

# From radvel but also found on StackExchange
class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())

lib_folder = os.path.dirname(os.path.realpath(__file__))
def get_requires():
    reqs = []
    for line in open('requirements.txt', 'r').readlines():
        reqs.append(line)
    return reqs

def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
                       open(project + '/__init__.py').read())
    return result.group(1)

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name = "ethraid",
    packages = find_packages(),
    include_package_data = True,
    package_data={'ethraid':['/data']},
    setup_requires=['cython', 'numpy'],
    install_requires=get_requires(),
    cmdclass={'build_ext':build_ext},
    entry_points = {
        "console_scripts": ['ethraid=ethraid.cli:main']
        },
    version = get_property('__version__', 'ethraid'),
    description = "Characterize long-period companions using RV trends, astrometric accelerations, and direct imaging",
    long_description = long_description,
    author = "Judah Van Zandt",
    author_email = "judahvz@astro.ucla.edu",
    url = "",
    ext_modules = cythonize(extensions, build_dir = 'ethraid/compiled', language_level='3', annotate=True)
    )