#ifndef VANDEVEN_H
#define VANDEVEN_H
#include "fourier_domain_filter.hpp"

class Vandeven: public FourierDomainFilter{
public:
    Vandeven(unsigned int order);
    virtual ~Vandeven();

    double evaluate(double omega) const;

    virtual void apply(MMSP::grid<2, MMSP::vector<fftw_complex> > &grid) const override;
    virtual void apply(MMSP::grid<3, MMSP::vector<fftw_complex> > &grid) const override;
private:
    unsigned int p{1};
    double *data{nullptr};

    /** Return x */
    inline double get_x(unsigned int indx) const;
    inline double dx() const;
    inline int index(double x) const;

    /** Initialise the data array */
    void precalculate_data();

    /** Evaluate the prefactor */
    double prefactor() const;

    template<int dim>
    void apply_generic(MMSP::grid<dim, MMSP::vector<fftw_complex> > &gr) const;
};

template<int dim>
void Vandeven::apply_generic(MMSP::grid<dim, MMSP::vector<fftw_complex> > &gr) const{
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