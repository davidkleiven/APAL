#include "test_build_matrices.hpp"
#include "tools.hpp"
#include "sparse_matrix.hpp"
#include "dense_matrix.hpp"
#include "use_numpy.hpp"

#include <cmath>

const unsigned int SIZE_3D_BUILD_MATRIX = 16;

PyObject* biharmonic_matrix(){
    SparseMatrix sp_mat;
    double alpha = 1.0;
    double M = 1.0;
    double dt = 1.0;

    cahn_hilliard_system_matrix3D(SIZE_3D_BUILD_MATRIX, M, alpha, dt, sp_mat);

    // Convert the sparse matrix to a full matrix
    unsigned int N = pow(SIZE_3D_BUILD_MATRIX, 3);
    DenseMatrix dense_matrix(N, N);

    sp_mat.to_dense(dense_matrix);

    // Convert dense matrix to numpy array
    npy_intp dims[2] = {N, N};
    PyObject* npy_array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    // Transfer dense-matrix to the Numpy array
    for (unsigned int row=0;row<N;row++)
    for (unsigned int col=0;col<N;col++){
        double* ptr = static_cast<double*>(PyArray_GETPTR2(npy_array, row, col));
        *ptr = dense_matrix(row, col);
    }
    return npy_array;
}