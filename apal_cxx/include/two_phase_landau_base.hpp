#ifndef TWO_PHASE_LANDAU_BASE_H
#define TWO_PHASE_LANDAU_BASE_H
#include<vector>

class TwoPhaseLandauBase{
public:
    TwoPhaseLandauBase(){};

    double evaluate(double conc, const std::vector<double> &shape) const;
    virtual double evaluate_vec(double x[]) const = 0;

    /** Evaluate the derivative */
    double partial_deriv_conc(double conc, const std::vector<double> &shape) const;
    virtual double partial_deriv_conc_vec(double x[]) const = 0;

    /** Partial derivative with respect to the shape variable */
    double partial_deriv_shape(double conc, const std::vector<double> &shape, unsigned int direction) const;
    virtual double partial_deriv_shape_vec(double x[], unsigned int direction) const = 0;
protected:
    static void conc_shape2vec(double conc, const std::vector<double> &shape, double vec[]);
};
#endif