#include "gaussian_filter.hpp"

void GaussianFilter::apply(MMSP::grid<2, MMSP::vector<fftw_complex> > &gr) const{
    apply_generic(gr);
}

void GaussianFilter::apply(MMSP::grid<3, MMSP::vector<fftw_complex> > &gr) const{
    apply_generic(gr);
}

double GaussianFilter::evaluate(double omega) const{
    // PI is defined in another file included in mmsp_files.cpp
    return exp(-pow(omega/width, 2));
}
