from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy
from distutils.sysconfig import get_python_lib
import sys

setup(name="lightcurve_cython",
      ext_modules=cythonize(Extension("lightcurve_cython",
                                      ["util.c", "lightcurve.pyx"],
                                      extra_compile_args=["-O3", "-fopenmp", "-std=c99"],
                                      extra_link_args=["-O3", "-fopenmp"],
                                      include_dirs=[numpy.get_include()]),
                            compiler_directives={'language_level': "3"},
                            annotate=True
                            )
      )
