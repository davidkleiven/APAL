#include "pfc1m.hpp"

using namespace std;
using cdouble = MMSP::scalar<complex<double> >;

template<int dim>
PFC1M<dim>::PFC1M(int L, const string &prefix, double eps, double alpha, double dt): PhaseFieldSimulation<dim>(L, prefix, 1), eps(eps), alpha(alpha), dt(dt){
    int dims[3] = {L, L, L};
    fft = new FFTW(dims, dim);

    if (dim == 1){
        cmplx_grid_ptr = new MMSP::grid<dim, cdouble>(0, L);
    } else if (dim == 2){
        cmplx_grid_ptr = new MMSP::grid<dim, cdouble>(0, L, 0, L);
    } else if (dim == 3){
        cmplx_grid_ptr = new MMSP::grid<dim, cdouble>(0, L, 0, L, 0, L);
    }
}

template<int dim>
PFC1M<dim>::~PFC1M(){
    delete fft;
    delete cmplx_grid_ptr;
}

template<int dim>
void PFC1M<dim>::to_parent_grid() {
    for (unsigned int i=0;i<MMSP::nodes(*this->grid_ptr);i++)
    {
        (*(this->grid_ptr))(i)[0] = real((*(this->cmplx_grid_ptr))(i));
        imag((*(this->cmplx_grid_ptr))(i)) = 0.0;
    }
}

template<int dim>
void PFC1M<dim>::from_parent_grid() {
    for (unsigned int i=0;i<MMSP::nodes(*this->grid_ptr);i++)
    {
        real((*(this->cmplx_grid_ptr))(i)) = (*(this->grid_ptr))(i)[0];
        imag((*(this->cmplx_grid_ptr))(i)) = 0.0;
    }
}

template<int dim>
void PFC1M<dim>::update(int nsteps) {
    #ifndef HAS_FFTW
        throw runtime_error("PFC1M requires FFTW!");
    #endif
    from_parent_grid();

    MMSP::grid<dim, cdouble>& gr = *(this->cmplx_grid_ptr);
    MMSP::grid<dim, cdouble> ft_cubed(gr);
    MMSP::grid<dim, cdouble> tmp_field(gr);

    for (unsigned int step=0;step<nsteps;step++){
        for (unsigned int node=0;node<MMSP::nodes(gr);node++) {
            real(tmp_field(node)) = pow(real(gr(node)), 3);
            imag(tmp_field(node)) = 0.0;
        }

        // Output into ft_cubed
        fft->execute(tmp_field, ft_cubed, FFTW_FORWARD);

        // Output fourier transform of fields into tmp_field
        fft->execute(gr, tmp_field, FFTW_FORWARD);

        for (unsigned int node=0;node<MMSP::nodes(gr);node++) {
            MMSP::vector<int> pos = gr.position(node);
            MMSP::vector<double> k_vec(pos.length());
            k_vector(pos, k_vec, this->L);
            double k = norm(k_vec);
            double k2 = pow(k, 2);
            double k4 = pow(k, 4);
            double dt = this->dt;

            real(tmp_field(node)) = real(tmp_field(node)) - dt*k2*real(ft_cubed(node))/(1.0 + dt*k2*(-eps + 1.0 + 2*k2 + k4));
            imag(tmp_field(node)) = imag(tmp_field(node)) - dt*k2*imag(ft_cubed(node))/(1.0 + dt*k2*(-eps + 1.0 + 2*k2 + k4));
        }
        fft->execute(tmp_field, *this->cmplx_grid_ptr, FFTW_BACKWARD);
    }
    to_parent_grid();
}

// Explicit instantiations
template class PFC1M<1>;
template class PFC1M<2>;
template class PFC1M<3>;