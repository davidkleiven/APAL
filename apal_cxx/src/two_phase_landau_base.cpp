#include "two_phase_landau_base.hpp"
#include <cstring>

using namespace std;

void TwoPhaseLandauBase::conc_shape2vec(double conc, const vector<double> &shape, double vec[]){
    vec[0] = conc;
    memcpy(vec+1, &shape[0], shape.size()*sizeof(double));
}

double TwoPhaseLandauBase::evaluate(double conc, const vector<double> &shape) const{
    double x[4];
    conc_shape2vec(conc, shape, x);
    return evaluate_vec(x);
}

double TwoPhaseLandauBase::partial_deriv_shape(double conc, const std::vector<double> &shape, unsigned int direction) const{
    double x[4];
    conc_shape2vec(conc, shape, x);
    return partial_deriv_shape_vec(x, direction);
}

double TwoPhaseLandauBase::partial_deriv_conc(double conc, const vector<double> &shape) const{
    double x[4];
    conc_shape2vec(conc, shape, x);
    return partial_deriv_conc_vec(x);
}