#include "dense_matrix.hpp"
#include <fstream>

using namespace std;

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

bool DenseMatrix::is_valid(unsigned int row, unsigned int col) const{
    return (row < nrows) && (col < ncols);
}

void DenseMatrix::save(const string &fname) const{
    ofstream out(fname);

    if (!out.good()){
        throw runtime_error("Could not open file for writing!");
    }

    for (unsigned int row=0;row<nrows;row++){
        for (unsigned int col=0;col<ncols;col++){
            out << (*this)(row, col);
            if (col < ncols-1){
                out << ",";
            }
        }
        out << "\n";
    }
}