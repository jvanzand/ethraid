from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize
from setuptools.command.build_ext import build_ext as _build_ext
import re

import os
import sys
print('FILE PATH', os.path.abspath(os.path.dirname(__file__)))
# import line_profiler
# profile = line_profiler.LineProfiler()

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
## include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]
extensions = [
    Extension("ethraid.compiled._kepler",
             ['ethraid/c_kepler/_kepler.pyx']),
    Extension("ethraid.compiled.helper_functions_general",
             ['ethraid/helper_functions_general.pyx']),
    Extension("ethraid.compiled.helper_functions_rv",
             ['ethraid/helper_functions_rv.pyx']),
    Extension("ethraid.compiled.helper_functions_astro",
             ['ethraid/helper_functions_astro.pyx'])
]

##############################

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
    setup_requires=['numpy'],
    install_requires=get_requires(),
    cmdclass={'build_ext':build_ext},
    entry_points = {
        "console_scripts": ['ethraid=ethraid.cli:main']
        },
    version = get_property('__version__', 'ethraid'),
    description = "Command line application for trend code",
    long_description = long_description,
    author = "Judah Van Zandt",
    author_email = "judahvz@astro.ucla.edu",
    url = "",
    ext_modules = cythonize(extensions, build_dir = 'ethraid/compiled', language_level='3', annotate=True)
    )

# ext_modules = cythonize(extensions, build_dir = 'ethraid/compiled', language_level='3', annotate=True)