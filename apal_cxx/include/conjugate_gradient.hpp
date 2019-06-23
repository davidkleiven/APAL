#ifndef CONJUGATE_GRADIENT_H
#define CONJUGATE_GRADIENT_H

#include "sparse_matrix.hpp"
#include <vector>

class ConjugateGradient{
public:
    ConjugateGradient(double tol);

    void solve(const SparseMatrix &matrix, const std::vector<double> &rhs, std::vector<double> &x) const;

    unsigned int get_num_iter() const {return last_num_iter;};
private:
    double tol{1E-6};
    unsigned int last_num_iter{0};
};
#endif