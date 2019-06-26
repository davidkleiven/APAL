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

PyObject* laplacian_matrix3D(){
    SparseMatrix sp_mat;
    double factor = 1.0;

    system_matrix_implicit_laplacian3D(SIZE_3D_BUILD_MATRIX, factor, sp_mat);

    unsigned int N = pow(SIZE_3D_BUILD_MATRIX, 3);
    DenseMatrix dense_matrix(N, N);
    sp_mat.to_dense(dense_matrix);

    // Convert to numpy object

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

PyObject* test_add_matrices_same_elements(){
    SparseMatrix mat1;
    SparseMatrix mat2;

    unsigned int rows[3] = {1, 4, 6};
    unsigned int cols[3] = {0, 2, 4};
    double values1[3] = {1.0, 2.0, 3.0};
    double values2[3] = {-2.0, 1.0, 2.3};

    unsigned int N = 7;
    DenseMatrix expected_sum(N, N);
    expected_sum.fill_zeros();
    for (unsigned int i=0;i<3;i++){
        mat1.insert(rows[i], cols[i], values1[i]);
        mat2.insert(rows[i], cols[i], values2[i]);
        expected_sum(rows[i], cols[i]) = values1[i] + values2[i];
    }

    mat1+= mat2;

    DenseMatrix summed(N, N);
    summed.fill_zeros();
    mat1.to_dense(summed);
    npy_intp dims[2] = {N, N};
    PyObject* npy_array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    PyObject* npy_array_expected = PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    // Transfer dense-matrix to the Numpy array
    for (unsigned int row=0;row<N;row++)
    for (unsigned int col=0;col<N;col++){
        double* ptr = static_cast<double*>(PyArray_GETPTR2(npy_array, row, col));
        *ptr = summed(row, col);

        ptr = static_cast<double*>(PyArray_GETPTR2(npy_array_expected, row, col));
        *ptr = expected_sum(row, col);
    }

    PyObject* npy_arrays = PyList_New(2);
    PyList_SetItem(npy_arrays, 0, npy_array);
    PyList_SetItem(npy_arrays, 1, npy_array_expected);
    return npy_arrays;
}

PyObject* test_add_matrices_mixed_elements(){
    SparseMatrix mat1;
    SparseMatrix mat2;

    unsigned int rows[3] = {1, 4, 6};
    unsigned int cols[3] = {0, 2, 4};
    unsigned int rows1[3] = {1, 2, 5};
    unsigned int cols1[3] = {0, 2, 3};
    double values1[3] = {1.0, 2.0, 3.0};
    double values2[3] = {-2.0, 1.0, 2.3};

    unsigned int N = 7;
    DenseMatrix expected_sum(N, N);
    expected_sum.fill_zeros();
    for (unsigned int i=0;i<3;i++){
        mat1.insert(rows[i], cols[i], values1[i]);
        mat2.insert(rows1[i], cols1[i], values2[i]);
        expected_sum(rows[i], cols[i]) += values1[i];
        expected_sum(rows1[i], cols1[i]) += + values2[i];
    }

    mat1+= mat2;

    DenseMatrix summed(N, N);
    summed.fill_zeros();
    mat1.to_dense(summed);
    npy_intp dims[2] = {N, N};
    PyObject* npy_array = PyArray_SimpleNew(2, dims, NPY_DOUBLE);
    PyObject* npy_array_expected = PyArray_SimpleNew(2, dims, NPY_DOUBLE);

    // Transfer dense-matrix to the Numpy array
    for (unsigned int row=0;row<N;row++)
    for (unsigned int col=0;col<N;col++){
        double* ptr = static_cast<double*>(PyArray_GETPTR2(npy_array, row, col));
        *ptr = summed(row, col);

        ptr = static_cast<double*>(PyArray_GETPTR2(npy_array_expected, row, col));
        *ptr = expected_sum(row, col);
    }

    PyObject* npy_arrays = PyList_New(2);
    PyList_SetItem(npy_arrays, 0, npy_array);
    PyList_SetItem(npy_arrays, 1, npy_array_expected);
    return npy_arrays;
}