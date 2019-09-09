#ifndef GAUSSIAN_FILTER_H
#define GAUSSIAN_FILTER_H
#include "fourier_domain_filter.hpp"

class GaussianFilter: public FourierDomainFilter{
public:
    GaussianFilter(double width): width(width){};

    virtual void apply(MMSP::grid<2, MMSP::vector<fftw_complex> > &grid) const override;
    virtual void apply(MMSP::grid<3, MMSP::vector<fftw_complex> > &grid) const override;

    double evaluate(double omega) const;
private:
    double width{0.0};

    template<int dim>
    void apply_generic(MMSP::grid<dim, MMSP::vector<fftw_complex> > &gr) const;
};


template<int dim>
void GaussianFilter::apply_generic(MMSP::grid<dim, MMSP::vector<fftw_complex> > &gr) const{
    double L = MMSP::xlength(gr);
    for (unsigned int node=0;node<MMSP::nodes(gr);node++){
        MMSP::vector<int> pos = gr.position(node);
        MMSP::vector<double> k_vec(pos.length());
        k_vector(pos, k_vec, L);

        double k = norm(k_vec);

        double filter_value = evaluate(k);

        // Multiply all fields with the filter value
        for (unsigned int field=0;field<MMSP::fields(gr);field++){
            real(gr(node)[field]) *= filter_value;
            imag(gr(node)[field]) *= filter_value;
        }
    }
};

#endif