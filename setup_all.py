from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext
from Cython.Build import cythonize
import re
import numpy as np
import sys



####################

# Extra stuff to profile code using kernprof -l -v file.py after it is compiled
# Taken from https://stackoverflow.com/questions/28301931/how-to-profile-cython-functions-line-by-line

from Cython.Compiler.Options import get_directive_defaults
directive_defaults = get_directive_defaults()

directive_defaults['linetrace'] = True
directive_defaults['binding'] = True

extensions = [
    Extension("helper_functions_wrapper", ['helper_functions_wrapper.pyx'], 
                include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')]),
    Extension("helper_functions", ['helper_functions.pyx'], 
                include_dirs=[np.get_include()], define_macros=[('CYTHON_TRACE', '1')])
]
##############################


class build_ext(_build_ext):
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process:
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


def get_property(prop, project):
    result = re.search(r'{}\s*=\s*[\'"]([^\'"]*)[\'"]'.format(prop),
                       open(project + '/__init__.py').read())
    return result.group(1)

path_to_setup = sys.path[0]
path_to_setup_dots = sys.path[0].strip('/').replace('/', '.')

path_list = [('c_kepler._kepler', path_to_setup+'/c_kepler/_kepler.pyx'),
             ('no_classes.helper_functions', path_to_setup+'/no_classes/helper_functions.pyx'),
             ('no_classes.helper_functions_wrapper', path_to_setup+'/no_classes/helper_functions_wrapper.pyx')]

print(path_list[0])


extensions = []
for i in path_list:
    extensions.append(eval('Extension("{}", ["{}"])'.format(i[0], i[1])))


#extensions = [Extension("_kepler", ["_kepler.pyx"],)]

reqs = []
for line in open(sys.path[0]+'/requirements.txt', 'r').readlines():

    reqs.append(line)

setup(
    packages=find_packages(),
    setup_requires=['numpy', 'cython'],
    ext_modules=cythonize(extensions, language_level=3, annotate=False),
    cmdclass={'build_ext': build_ext},
    install_requires=reqs,
    include_package_data=True
)