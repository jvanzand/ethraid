from setuptools import setup
from Cython.Build import cythonize
from distutils.extension import Extension

# Extension(name_of_so_file, [name_of_pyx_file.pyx])
# Because of the trends/ at the beginning of the path, this needs to be run from the top-level directory.
# Is there a way to allow it to be run from both? Maybe, check with Erik.
extensions = [
    Extension("example_cython", ["example_cython.pyx"])
]

setup(
    ext_modules = cythonize(extensions, language_level='3', annotate=True)
)