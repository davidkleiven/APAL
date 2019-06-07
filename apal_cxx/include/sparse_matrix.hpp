#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H
#include <vector>
#include <string>

// Forward-declaration
class DenseMatrix;

class SparseMatrix{
public:
    SparseMatrix(){};

    void clear();
    void insert(unsigned int row, unsigned int col, double value);

    void dot(const std::vector<double> &vec, std::vector<double> &res) const;

    void save(const std::string &fname) const;

    bool is_symmetric() const;

    /** Count the number of entries on each row */
    void count_entries_on_each_row(std::vector<unsigned int > &num_entries) const;

    /** Convert to CSR format */
    void to_csr();

    /** Convert to dense matrix. It is assumed that dense_mat has the correct number of rows and columns */
    void to_dense(DenseMatrix &dense_mat) const;

    /** Check if the same row column pair is inserted multiple times */
    bool has_duplicates() const;
private:
    bool converted_to_csr{false};
    unsigned int num_rows{0};
    std::vector<double> values;
    std::vector<unsigned int> row;
    std::vector<unsigned int> col;

    // CSR related structures
    std::vector<double> csr_values;
    std::vector<int> csr_col_indx;
    std::vector<int> row_ptr;

    /** Return the largest row index */
    unsigned int max_row() const;
};
#endif