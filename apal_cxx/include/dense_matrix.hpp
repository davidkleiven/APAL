#ifndef DENSE_MATRIX_H
#define DENSE_MATRIX_H
#include <string>

class DenseMatrix{
public:
    DenseMatrix(unsigned int nrows, unsigned int ncols);
    ~DenseMatrix();

    /** Read write access operator */
    double& operator()(int row, int col){return data[index(row, col)];};

    /** Read-only access operator */
    const double& operator()(int row, int col) const {return data[index(row, col)];};

    /** Fill with zeros */
    void fill_zeros();

    /** Check if row,col pair is valid */
    bool is_valid(unsigned int row, unsigned int col) const;

    /** Store matrix in CSV format */
    void save(const std::string &fname) const;
private:
    unsigned int nrows{1};
    unsigned int ncols{1};
    unsigned int index(unsigned int row, unsigned int col) const;
    double *data{nullptr};
};
#endif