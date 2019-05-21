import sys
from setuptools import setup, Extension, find_packages
import distutils
from distutils.errors import CompileError, LinkError
import numpy as np
from Cython.Build import cythonize
import os
from textwrap import dedent

fftw_libs = ["fftw3_threads", "fftw3"]

def check_fftw():
    """Check if FFTW extension works."""

    cpp_code = dedent("""
        #include <fftw.h>

        int main(int argc, char *argv[]){
            return 0;
        }
    """)

    # Store the code to a temporary file
    fname = "check_fftw.cpp"
    binfname = "check_fftw.out"
    with open(fname, 'w') as out:
        out.write(cpp_code)
    
    compiler = distutils.ccompiler.new_compiler()
    ret_val = True
    try:
        compiler.link_executable(
            compiler.compile([fname]),
            binfname,
            libraries=fftw_libs,
        )
    except CompileError:
        print('FFTW compile error')
        ret_val = False
    except LinkError:
        print('FFTW link error')
        ret_val = False
    
    os.remove(fname)
    os.remove(binfname)
    os.remove(fname.split(".")[0] + ".o")
    return ret_val


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

if check_fftw():
    optional_lib_phasefield += fftw_libs
    define_macros.append(("HAS_FFTW", 1))
    print("Compiling with FFTW")

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