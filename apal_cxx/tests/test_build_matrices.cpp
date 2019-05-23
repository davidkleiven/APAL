#include "test_build_matrices.hpp"
#include "tools.hpp"
#include "sparse_matrix.hpp"
#include "dense_matrix.hpp"

#include <cmath>

const unsigned int L = 64;

PyObject* biharmonic_matrix(){
    SparseMatrix sp_mat;
    double alpha = 1.0;
    double M = 1.0;
    double dt = 1.0;

    cahn_hilliard_system_matrix3D(L, M, alpha, dt, mat);

    // Convert the sparse matrix to a full matrix
    unsigned int N = pow(L, 3);
    DenseMatrix dense_matrix(N, N);

    mat.to_dense(dense_matrix);

    // Convert dense matrix to numpy array
    npy_intp dims[2] = {N, N};
    PyObject* npy_array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    // Transfer dense-matrix to the Numpy array
    for (unsigned int row=0;row<N;row++)
    for (unsigned int col=0;col<N;col++){
        double* ptr = static_cast<double*>(PyArray_GETPTR2(npy_array, row, col));
        *ptr = mat(row, col);
    }
    return npy_array;
}