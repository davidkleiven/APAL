#include "test_build_matrices.hpp"
#include "tools.hpp"
#include "sparse_matrix.hpp"

const unsigned int L = 64;

PyObject* biharmonic_matrix(){
    SparseMatrix sp_mat;
    double alpha = 1.0;
    double M = 1.0;
    double dt = 1.0;

    cahn_hilliard_system_matrix3D(L, M, alpha, dt, mat);

    // Convert the sparse matrix to a full matrix

}