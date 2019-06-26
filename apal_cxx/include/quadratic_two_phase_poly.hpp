#ifndef QUADRATIC_TWO_PHASE_H
#define QUADRATIC_TWO_PHASE_H

#include "two_phase_landau_base.hpp"
#include "polynomial.hpp"

class QuadraticTwoPhase: public TwoPhaseLandauBase{
public:
    QuadraticTwoPhase(){};

    /** Set polynomial for the first phase */
    void set_poly_phase1(const Polynomial &poly);

    /** Set polynomial for the second phase */
    void set_poly_phase2(const Polynomial &poly);

    virtual void in_valid_state() const override;

    virtual double evaluate_vec(double x[]) const override;

    virtual double partial_deriv_conc_vec(double x[]) const override;

    virtual double partial_deriv_shape_vec(double x[], unsigned int direction) const override;
private:
    const Polynomial *phase1{nullptr};
    const Polynomial *phase2{nullptr};
};
#endif