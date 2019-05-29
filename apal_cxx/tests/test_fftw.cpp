#include "test_fftw.hpp"
#include <Python.h>
#include "fftw_mmsp.hpp"
#include "fftw_complex_placeholder.hpp"
#include <cmath>
#include <iostream>

using namespace std;

void fft1D(vector<double> &array, int direction){
    #ifndef HAS_FFTW
        return Py_NONE;
    #endif

    int dims[1] = {array.size()};
    FFTW fftw(1, dims);
    vector<fftw_complex> out(array.size());
    fftw.execute(array, out, direction);

    for (unsigned int i=0;i<out.size();i++){
        array[i] = real(out[i]);
    }
}


void fft2D(vector<double> &array, int direction){
    #ifndef HAS_FFTW
        return Py_NONE;
    #endif

    int N = static_cast<int>(sqrt(array.size()));
    int dims[2] = {N, N};
    FFTW fftw(2, dims);
    vector<fftw_complex> out(array.size());
    fftw.execute(array, out, direction);

    for (unsigned int i=0;i<out.size();i++){
        array[i] = real(out[i]);
    }
}

void fft3D(vector<double> &array, int direction){
    #ifndef HAS_FFTW
        return Py_NONE;
    #endif

    int N = static_cast<int>(round(pow(array.size(), 1.0/3.0)));
    int dims[3] = {N, N, N};
    FFTW fftw(3, dims);
    vector<fftw_complex> out(array.size());
    fftw.execute(array, out, direction);

    for (unsigned int i=0;i<out.size();i++){
        array[i] = real(out[i]);
    }
}