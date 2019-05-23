#include "dense_matrix.hpp"

DenseMatrix::DenseMatrix(unsigned int nrows, unsigned int ncols): nrows(nrows), ncols(ncols){
    data = new double[nrows*ncols];
}

DenseMatrix::~DenseMatrix(){
    delete [] data;
}

unsigned int DenseMatrix::index(unsigned int row, unsigned int col) const{
    return row*ncols + col;
}

void DenseMatrix::fill_zeros(){
    for (unsigned int i=0;i<nrows*ncols;i++){
        data[i] = 0.0;
    }
}