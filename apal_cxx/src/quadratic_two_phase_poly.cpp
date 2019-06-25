#include "quadratic_two_phase_poly.hpp"
#include <stdexcept>

using namespace std;

void QuadraticTwoPhase::in_valid_state() const{
    if (!phase1){
        throw invalid_argument("Polynomial for phase1 is not set!");
    }

    if (!phase2){
        throw invalid_argument("Polynomial for phase2 is not set!");
    }
}

void QuadraticTwoPhase::set_poly_phase1(const Polynomial &poly){
    if (!poly.get_dim() != 1){
        throw invalid_argument("Polynomial for the first phase can only depend on one variable!");
    }
    phase1 = &poly;
}

void QuadraticTwoPhase::set_poly_phase2(const Polynomial &poly){
    phase2 = &poly;
}

double QuadraticTwoPhase::evaluate_vec(double x[]) const{
    return phase1->evaluate(x) + phase2->evaluate(x);
}

double QuadraticTwoPhase::partial_deriv_conc_vec(double x[]) const{
    return phase1->deriv(x, 0) + phase2->deriv(x, 0);
}

double QuadraticTwoPhase::partial_deriv_shape_vec(double x[], unsigned int direction) const{
    return phase2->deriv(x, direction+1);
}


