import sys
from setuptools import setup, Extension, find_packages
import numpy as np
from Cython.Build import cythonize
import os


src_phase = "apal_cxx/src"
phasefield_sources = ["mat4D.cpp", "khacaturyan.cpp",
                      "linalg.cpp", "cahn_hilliard.cpp", "polynomial_term.cpp",
                      "polynomial.cpp", "regression_kernels.cpp", "kernel_regressor.cpp",
                      "sparse_matrix.cpp", "gaussian_white_noise.cpp",
                      "adaptive_timestep_logger.cpp", "additional_tools.cpp"]

phasefield_sources = [src_phase + "/" + x for x in phasefield_sources]
phasefield_sources.append("apal/cython/apal_cxx.pyx")

extra_comp_args = []
define_macros = []
optional_lib_phasefield = []
phase_field_mod = Extension("apal_cxx", sources=phasefield_sources,
                            include_dirs=[np.get_include(),
                                          "apal_cxx/include",
                                          "apal_cxx/src",
                                          "apal_cxx/tests",
                                          os.environ.get("MMSP_HOME", "./")+"/include"],
                            extra_compile_args=extra_comp_args,
                            language="c++", define_macros=define_macros,
                            libraries=["gomp", "pthread", "z"] + optional_lib_phasefield)

setup(
    name="apal",
    ext_modules=cythonize(phase_field_mod),
    version=1.0,
    author="David Kleiven",
    author_email="davidkleiven446@gmail.com",
    description="Alloy Phasefield Abstraction Layer",
    packages=find_packages(),
    include_package_data=True
)