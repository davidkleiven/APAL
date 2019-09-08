#ifndef FOURIER_DOMAIN_FILTER_H
#define FOURIER_DOMAIN_FILTER_H

#include "MMSP.grid.h"
#include "MMSP.vector.h"
#include "fftw_complex_placeholder.hpp"

class FourierDomainFilter{
public:
    virtual ~FourierDomainFilter(){};

    /** Apply method */
    virtual void apply(MMSP::grid<2, MMSP::vector<fftw_complex> > &grid) const = 0;
    virtual void apply(MMSP::grid<3, MMSP::vector<fftw_complex> > &grid) const = 0;
};
#endif