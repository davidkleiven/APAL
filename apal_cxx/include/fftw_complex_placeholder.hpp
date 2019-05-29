#ifndef FFTW_COMPLEX_PLACEHOLDER_H
#define FFTW_COMPLEX_PLACEHOLDER_H

#ifndef HAS_FFTW
#include <complex>

struct apal_complex_t{
    double re;
    double im;
};

//typedef double fftw_complex[2];
typedef apal_complex_t fftw_complex;
typedef int fftw_direction;

const int FFTW_FORWARD = 1;
const int FFTW_BACKWARD = -1;
const int FFTW_ESTIMATE = 0x01;

// Dummy value for an fftw ndplan
struct fftwnd_plan{};
struct fftw_plan{};
void fftw_execute(fftw_plan &plan){};
#endif
#endif