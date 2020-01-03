from Cython.Build import cythonize
from distutils.core import setup
from distutils.extension import Extension
from numpy import get_include

ext_modules = cythonize(['pyqxmath.pyx'])
setup(ext_modules=ext_modules, include_dirs=[get_include()])