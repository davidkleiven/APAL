#include "tools.hpp"
#include "origin_singularity_integration.hpp"
#include <cmath>
#include <omp.h>
#include <stdexcept>

const double PI = acos(-1.0);

using namespace std;

void k_vector(const MMSP::vector<int> &pos, MMSP::vector<double> &k_vec, int N){
    for (unsigned int i=0;i<pos.length();i++){
        if (pos[i] < N/2){
            k_vec[i] = PI*pos[i]/N;
        }
        else{
            k_vec[i] = -PI*(N-pos[i])/N;
        }
    }
}

void dot(const mat3x3 &mat1, const MMSP::vector<double> &vec, MMSP::vector<double> &out){
    for (unsigned int i=0;i<3;i++){
        out[i] = 0.0;
    }

    for (unsigned int i=0;i<3;i++)
    for (unsigned int j=0;j<3;j++){
        out[i] += mat1[i][j]*vec[j];
    }
}

double dot(const MMSP::vector<double> &vec1, const MMSP::vector<double> &vec2){
    double value = 0.0;
    for (unsigned int i=0;i<vec1.length();i++){
        value += vec1[i]*vec2[i];
    }
    return value;
}

double dot(const vector<double> &v1, const vector<double> &v2){
    double inner_prod = 0.0;
    #ifndef NO_PHASEFIELD_PARALLEL
    #pragma omp parallel for reduction(+ : inner_prod)
    #endif
    for (unsigned int i=0;i<v1.size();i++){
        inner_prod += v1[i]*v2[i];
    }
    return inner_prod;
}

double norm(const MMSP::vector<double> &vec){
    double value = 0.0;
    for (unsigned int i=0;i<vec.length();i++){
        value += pow(vec[i], 2);
    }
    return sqrt(value);
}

void inplace_minus(vector<double> &vec1, const vector<double> &vec2){
    #ifndef NO_PHASEFIELD_PARALLEL
    #pragma omp parallel for
    #endif
    for (unsigned int i=0;i<vec1.size();i++){
        vec1[i] -= vec2[i];
    }
}

double inf_norm(const vector<double> &vec){
    double max_val = abs(vec[0]);
    #ifndef NO_PHASEFIELD_PARALLEL
    #pragma omp parallel for reduction(max : max_val)
    #endif
    for (unsigned int i=0;i<vec.size();i++){
        max_val = abs(vec[i]) > max_val ? abs(vec[i]) : max_val;
    }
    return max_val;
}

double least_squares_slope(double x[], double y[], unsigned int N){
    double Sxx = 0.0;
    double Sxy = 0.0;
    double Sx = 0.0;
    double Sy = 0.0;

    for (unsigned int i=0;i<N;i++){
        Sx += x[i];
        Sy += y[i];
        Sxy += x[i]*y[i];
        Sxx += x[i]*x[i];
    }

    return (N*Sxy - Sx*Sy)/(N*Sxx - Sx*Sx);
}

double contract_tensors(const mat3x3 &mat1, const mat3x3 &mat2){
    double value = 0.0;
    for (unsigned int i=0;i<3;i++)
    for (unsigned int j=0;j<3;j++){
        value += mat1[i][j]*mat2[i][j];
    }
    return value;
}

bool is_origin(MMSP::vector<double> &dir){
    double tol = 1E-6;
    for (int i=0;i<dir.length();i++){
        if (abs(dir[i]) > tol) return false;
    }
    return true;
}

double B_tensor_element(MMSP::vector<double> &dir, const mat3x3 &green, \
                        const mat3x3 &eff_stress1, const mat3x3 &eff_stress2)
{
    // if (is_origin(dir)){
    //     OriginSingularityIntegration integrator(10);
    //     vector< MMSP::vector<double> > directions;
    //     // TODO: Pass dimension as argument to support more than 2D
    //     integrator.get_integration_directions(2, directions);
    //     return B_tensor_element_origin(green, eff_stress1, eff_stress2, directions);
    // }

    MMSP::vector<double> temp_vec(3);
    dot(eff_stress2, dir, temp_vec);

    MMSP::vector<double> temp_vec2(3);
    dot(green, temp_vec, temp_vec2);
    dot(eff_stress1, temp_vec2, temp_vec);
    return dot(dir, temp_vec);
}

double B_tensor_element_origin(const mat3x3 &green, const mat3x3 &eff_stress1, const mat3x3 &eff_stress2, \
                               std::vector< MMSP::vector<double> > &directions)
{
    double value = 0.0;
    for (auto& dir : directions){
        value += B_tensor_element(dir, green, eff_stress1, eff_stress2);
    }
    return value/directions.size();
}

MMSP::vector<int>& wrap(MMSP::vector<int> &pos, unsigned int L){
    for (unsigned int i=0;i<pos.length();i++){
        if (pos[i] < 0){
            pos[i] += L;
        }
        else if (pos[i] >= L){
            pos[i] -= L;
        }
    }
    return pos;
}

void cahn_hilliard_system_matrix3D(unsigned int L, double M, double alpha, double dt, SparseMatrix &mat){
  mat.clear();

  // Create a grid containing the indices
  MMSP::grid<3, int> indexGrid(1, 0, L, 0, L, 0, L);
  for (unsigned int i=0;i<MMSP::nodes(indexGrid);i++){
      indexGrid(i) = i;
  }

  double factor = 2.0*alpha*M*dt;

  for (unsigned int i=0;i<MMSP::nodes(indexGrid);i++){

    // Insert diagonal
    mat.insert(i, i, 1.0 + 42*factor);

    // Insert elements that is one off
    // Retrive node at position +- 1
    for (unsigned int dir=0;dir<3;dir++)
    for (int j=-1;j<2;j+=2){
        MMSP::vector<int> pos = indexGrid.position(i);
        pos[dir] += j;
        unsigned int col = indexGrid(wrap(pos, L));
        mat.insert(i, col, -12*factor);
    }

    // Insert element 2 off the diagonal
    // Calculate factor at position +- 2
    for (unsigned int dir=0;dir<3;dir++)
    for (int j=-2;j<5;j+=4){
      MMSP::vector<int> pos = indexGrid.position(i);
      pos[dir] += j;
      unsigned int col = indexGrid(wrap(pos, L));
      mat.insert(i, col, factor);
    }    

    // Calculate the cross terms
    for (int ix=-1;ix<2;ix++)
    for (int iy=-1;iy<2;iy++)
    for (int iz=-1;iz<2;iz++)
    {
        // One of ix, iy, iz has to be zero
        if (abs(ix) + abs(iy) + abs(iz) != 2){
            continue;
        }

        MMSP::vector<int> pos = indexGrid.position(i);
        pos[0] += ix;
        pos[1] += iy;
        pos[2] += iz;
        unsigned int col = indexGrid(wrap(pos, L));
        mat.insert(i, col, 2*factor);
    }
  }

  if (!mat.is_symmetric()){
      throw runtime_error("System matrix 3D has to be symmetric!");
  }
  mat.to_csr();
}

void system_matrix_implicit_laplacian3D(unsigned int L, double prefactor[3], SparseMatrix &mat){
    mat.clear();

    // Create a grid containing the indices
    MMSP::grid<3, int> indexGrid(1, 0, L, 0, L, 0, L);
    for (unsigned int i=0;i<MMSP::nodes(indexGrid);i++){
        indexGrid(i) = i;
    }

    double sum_prefactor = prefactor[0] + prefactor[1] + prefactor[2];
    for (unsigned int i=0;i<MMSP::nodes(indexGrid);i++){
            // Insert diagonal
            mat.insert(i, i, 1.0 + 2*sum_prefactor);

            // Insert elements that is one off
            // Retrive node at position +- 1
            for (unsigned int dir=0;dir<3;dir++)
            for (int j=-1;j<2;j+=2){
                MMSP::vector<int> pos = indexGrid.position(i);
                pos[dir] += j;
                unsigned int col = indexGrid(wrap(pos, L));
                mat.insert(i, col, -prefactor[dir]);
            }
    }

    if (!mat.is_symmetric()){
        throw runtime_error("Laplacian 3D matrix has to be symmetric!");
    }
    mat.to_csr();
}

void system_matrix_implicit_laplacian3D(unsigned int L, double prefactor, SparseMatrix &mat){
    double pref_array[3] = {prefactor, prefactor, prefactor};
    system_matrix_implicit_laplacian3D(L, pref_array, mat);
}

bool is_origin(const MMSP::vector<double> &vec){
    const double tol = 1E-8;
    for (unsigned int i=0;i<vec.length();i++){
        if (abs(vec[i]) > tol){
            return false;
        }
    }
    return true;
}