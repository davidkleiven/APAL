#include "raised_cosine.hpp"
#include <cmath>

using namespace std;

double RaisedCosine::evaluate(double omega) const{
    if (omega < omega_min()){
        return 1.0;
    }
    else if (omega >= omega_max()){
        return 0.0;
    }

    // PI is defined in another file included in mmsp_files.cpp
    return 0.5*(1.0 + cos(PI*(omega - omega_min())/(2*roll_off*omega_cut)));
}

double RaisedCosine::omega_min() const{
    return (1.0 - roll_off)*omega_cut;
}

double RaisedCosine::omega_max() const{
    return (1.0 + roll_off)*omega_cut;
}

void RaisedCosine::apply(MMSP::grid<2, MMSP::vector<fftw_complex> > &gr) const{
    apply_generic(gr);
}

void RaisedCosine::apply(MMSP::grid<3, MMSP::vector<fftw_complex> > &gr) const{
    apply_generic(gr);
}