from distutils.core import setup
from distutils.extension import Extension
import numpy as np
from Cython.Distutils import build_ext


setup(
    name = 'Cython functions and classes',
    ext_modules = [Extension("cyFuns", ["cyFuns.pyx"])],
    cmdclass = {'build_ext': build_ext},
    include_dirs = [np.get_include()]
)
