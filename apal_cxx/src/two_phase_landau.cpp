#include "two_phase_landau.hpp"
#include "tools.hpp"
#include <cmath>

using namespace std;
TwoPhaseLandau::TwoPhaseLandau(): TwoPhaseLandauBase(){};

double TwoPhaseLandau::evaluate_vec(double x[]) const{
    double value = regressor->evaluate(x[0]);
    value += polynomial->evaluate(x);
    return value + jump_contrib(x[0]);
}


double TwoPhaseLandau::partial_deriv_conc_vec(double x[]) const{
    const unsigned int N = 5;
    const double conc_width = 0.05;
    const double dx = conc_width/(N-1);
    double current_conc = x[0]-conc_width/2.0;
    double concs[N];
    double energies[N];

    for (unsigned int i=0;i<N;i++){
        concs[i] = current_conc;
        x[0] = current_conc;
        current_conc += dx;
        energies[i] = this->evaluate_vec(x);
    }
    return least_squares_slope(concs, energies, N);


    double value = regressor->deriv(x[0]);
    value += polynomial->deriv(x, 0); 
    return value + jump_contrib_deriv(x[0]);
}

double TwoPhaseLandau::partial_deriv_shape_vec(double x[], unsigned int direction) const{
    return polynomial->deriv(x, direction+1);
}

unsigned int TwoPhaseLandau::get_poly_dim() const{
    if (!polynomial){
        throw runtime_error("Dimension of polynmoial requested, but no polynomial is set!");
    }
    return polynomial->get_dim();
}

void TwoPhaseLandau::set_discontinuity(double conc, double jump){
    disc_conc = conc;
    disc_jump = jump;
}

double TwoPhaseLandau::jump_contrib(double conc) const{
    double x = (conc - disc_conc)/step_func_width;
    return 0.0;
    return -disc_jump*0.5*(1 + tanh(x));
}

double TwoPhaseLandau::jump_contrib_deriv(double conc) const{
    double x = (conc - disc_conc)/step_func_width;
    return 0.0;
    return -disc_jump*0.5*(1 - pow(tanh(x), 2))/step_func_width;
}

void TwoPhaseLandau::in_valid_state() const{
    if (!get_regressor()){
        throw invalid_argument("TwoPhaseLanday has no kernel regressor!");
    }

    if (!get_regressor()->kernel_is_set()){
        throw invalid_argument("The Kernel Regressor has no kernel!");
    }

    if (!polynomial){
        throw invalid_argument("Polyomial is not set!");
    }
}