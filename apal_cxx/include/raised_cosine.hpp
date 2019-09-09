#ifndef RAISED_COSINE_H
#define RAISED_COSINE_H
#include "fourier_domain_filter.hpp"

class RaisedCosine: public FourierDomainFilter{
public:
    RaisedCosine(double omega_cut, double roll_off): omega_cut(omega_cut), roll_off(roll_off){};

    /** Evaluate the filter function */
    double evaluate(double omega) const;

    /** Return minimum angular frequency before filtering kicks in */
    double omega_min() const;

    /** Return the frequency beyond which all frequencies are damped */
    double omega_max() const;

    virtual void apply(MMSP::grid<2, MMSP::vector<fftw_complex> > &grid) const override;
    virtual void apply(MMSP::grid<3, MMSP::vector<fftw_complex> > &grid) const override;
private:
    double omega_cut{0.0};
    double roll_off{0.0};

    template<int dim>
    void apply_generic(MMSP::grid<dim, MMSP::vector<fftw_complex> > &gr) const;
};

template<int dim>
void RaisedCosine::apply_generic(MMSP::grid<dim, MMSP::vector<fftw_complex> > &gr) const{
    double L = MMSP::xlength(gr);
    for (unsigned int node=0;node<MMSP::nodes(gr);node++){
        MMSP::vector<int> pos = gr.position(node);
        MMSP::vector<double> k_vec(pos.length());
        k_vector(pos, k_vec, L);

        double filter_value = 1.0;
        for (unsigned int comp=0;comp<dim;comp++){
            filter_value *= evaluate(k_vec[comp]);
        }

        // Multiply all fields with the filter value
        for (unsigned int field=0;field<MMSP::fields(gr);field++){
            real(gr(node)[field]) *= filter_value;
            imag(gr(node)[field]) *= filter_value;
        }
    }
};
#endif