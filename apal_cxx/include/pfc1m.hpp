#ifndef PFC1M_H
#define PFC1M_H
#include <complex>
#include "phase_field_simulation.hpp"
#ifdef HAS_FFTW
    #include <fftw3.h>
#endif
#include "fftw_complex_placeholder.hpp"

template<int dim>
class PFC1M: public PhaseFieldSimulation<dim> {
    public:
        PFC1M(int L, const std::string &prefix, double eps, double alpha, double dt);
        virtual ~PFC1M();

        // Propagate one time step
        virtual void update(int nsteps) override;
    protected:
        double eps{1.0};
        double alpha{1.0};
        double dt{1.0};
        FFTW *fft{nullptr}

        MMSP::grid<dim, MMSP::scalar<std::complex<double> > > *cmplx_grid_ptr{nullptr};
        void to_parent_grid();
        void from_parent_grid();
}
#endif