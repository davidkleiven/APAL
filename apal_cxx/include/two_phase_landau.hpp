#ifndef TWO_PHASE_LANDAU_H
#define TWO_PHASE_LANDAU_H
#include <vector>
#include "kernel_regressor.hpp"
#include "polynomial.hpp"
#include "two_phase_landau_base.hpp"

class TwoPhaseLandau: public TwoPhaseLandauBase{
    public:
        TwoPhaseLandau();

        /** Set the kernel regressor */
        void set_kernel_regressor(const KernelRegressor &regr){regressor = &regr;};

        /** Set the polymial */
        void set_polynomial(const Polynomial &poly){polynomial = &poly;};

        /** Set discontinuity jump */
        void set_discontinuity(double conc, double jump);

        /** Get the dimension of the polynomial */
        unsigned int get_poly_dim() const;

        /** Evaluate the shape polynomial*/
        double evaluate_vec(double x[]) const override;

        /** Evaluate the derivative */
        double partial_deriv_conc_vec(double x[]) const override;

        /** Partial derivative with respect to the shape variable */
        double partial_deriv_shape_vec(double x[], unsigned int direction) const override;

        /** Return a pointer to the regressor */
        const KernelRegressor* get_regressor() const{return regressor;};
    private:
        const KernelRegressor *regressor{nullptr};
        const Polynomial *polynomial{nullptr};

        /** Return contribution from discontiuity jump */
        //double jump_contrib(double conc) const {return conc >= disc_conc ? -disc_jump : 0.0;};
        double jump_contrib(double conc) const;
        double jump_contrib_deriv(double conc) const;

        double disc_jump{0.0};
        double disc_conc{0.0};
        double step_func_width{0.1};
};
#endif