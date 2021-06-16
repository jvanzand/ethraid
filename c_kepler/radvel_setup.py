from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext as _build_ext
import re

"""
Modified by Judah to make changes to Kepler solver
"""

from Cython.Build import cythonize


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


extensions = [Extension("_kepler", ["_kepler.pyx"],)]

reqs = []
for line in open('radvel_requirements.txt', 'r').readlines():
    # if not line.startswith('celerite') and not line.startswith('h5py'):
    if not line.startswith('h5py'):
        reqs.append(line)

setup(
    packages=find_packages(),
    setup_requires=['numpy', 'cython'],
    ext_modules=cythonize(extensions, language_level=3, annotate=False),
    cmdclass={'build_ext': build_ext},
    data_files=[
        (
            'radvel_example_data', 
            [
                'example_data/164922_fixed.txt', 
                'example_data/epic203771098.csv',
                'example_data/k2-131.txt'
            ]
        )
    ],
    entry_points={'console_scripts': ['radvel=radvel.cli:main']},
    install_requires=reqs,
    include_package_data=True
)
