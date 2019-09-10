#include "vandeven.hpp"
#include "tools.hpp"

using namespace std;

const unsigned int NUM_TABULATED_VALUES = 1000;

Vandeven::Vandeven(unsigned int order): p(order){
    data = new double[NUM_TABULATED_VALUES];
    precalculate_data();
}

Vandeven::~Vandeven(){
    delete [] data;
}

void Vandeven::precalculate_data(){
    double pref = prefactor();

    // Calculate the integral
    data[0] = 0.0;
    for (unsigned int i=1;i<NUM_TABULATED_VALUES;i++){
        double x1 = get_x(i-1);
        double x2 = get_x(i);

        double y1 = pow(x1*(1-x1), p-1);
        double y2 = pow(x2*(1-x2), p-1);

        data[i] = data[i-1] + 0.5*(y1 + y2)*dx();
    }

    // Update the data array to actually hold the filter values
    for (unsigned int i=0;i<NUM_TABULATED_VALUES;i++){
        data[i] = 1.0 - pref*data[i];
    }
}

double Vandeven::dx() const{
    return 1.0/(NUM_TABULATED_VALUES-1);
}

double Vandeven::get_x(unsigned int i) const{
    return i*dx();
}

int Vandeven::index(double x) const{
    return x/dx();
}

double Vandeven::prefactor() const{
    // Calculate 
    double ln1 = ln_factorial(2*p-1);
    double ln2 = ln_factorial(p-1);
    double ratio = ln1 - 2*ln2;
    return exp(ratio);
}

double Vandeven::evaluate(double omega) const{
    double omega_max = PI/2.0;
    double scaled_omega = omega/omega_max;
    int i = index(scaled_omega);

    double y1 = data[i];
    double y2 = data[i];
    if (i < NUM_TABULATED_VALUES-1){
        y2 = data[i+1];
    }

    double w = (scaled_omega - get_x(i))/dx();
    return (1-w)*y1 + w*y2;
}

void Vandeven::apply(MMSP::grid<2, MMSP::vector<fftw_complex> > &gr) const{
    apply_generic(gr);
}

void Vandeven::apply(MMSP::grid<3, MMSP::vector<fftw_complex> > &gr) const{
    apply_generic(gr);
}