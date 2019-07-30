#ifndef CHGL_REALSPACE_H
#define CHGL_REALSPACE_H
#include "chgl.hpp"
#include "sparse_matrix.hpp"
#include <array>
#include <set>

template<int dim>
class CHGLRealSpace: public CHGL<dim>{
public:
    CHGLRealSpace(int L, const std::string &prefix, unsigned int num_gl_fields, \
         double M, double alpha, double dt, double gl_damping, 
         const interface_vec_t &interface);

    virtual ~CHGLRealSpace();

    /** Set how often the field derivatives should be updated */
    void set_field_update_rate(unsigned int rate);

    /** Set how often the strain fields should be updated */
    void set_strain_update_rate(unsigned int rate);

    /** Build the matrix for CHGL equations */
    void build2D();

    /** Build the matrix for CHGL equations in 3D */
    void build3D();

    /** Implement the update function */
    virtual void update(int nsteps) override;

    /** Calculate the energy of the system */
    void energy(std::map<std::string, double> &tr_item) const;

    /** Conserve volume of the order parameter squared */
    void conserve_volume(unsigned int gl_field);
private:
    unsigned int implicitDir{0};
    unsigned int field_deriv_update_freq{1};
    unsigned int strain_deriv_update_freq{1};
    std::array<SparseMatrix, dim+1> matrices;
    bool did_build_matrices{false};

    double min_strain_deriv{0.0};
    double max_strain_deriv{0.0};

    MMSP::vector<int> & wrap(MMSP::vector<int> &pos) const;
    std::set<unsigned int> conserved_gl_fields;

    unsigned int node_index(MMSP::vector<int> &pos) const;

    void add_cook_noise_to_fd_scheme(std::vector<double> &rhs, int field) const;
    void add_strain_contribution(std::vector<double> &rhs, int field) const;
    void calculate_strain_contribution();
    void update_min_max_strain_deriv();
    MMSP::grid<dim, MMSP::vector<fftw_complex> > *strain_deriv{nullptr};

    /** Return the prefactor for the central stencil of the laplacian */
    double laplacian_central_stencil() const;

    /** Check if the timestep should be lowered */
    bool should_lower_timestep(double energy) const;

    /** Log a track item */
    void log_tritem(const std::map<std::string, double> &item) const;

    /** Rebuild the system matrices */
    void rebuild_matrices();

    /** Write the strain derivative to file */
    void save_strain_derivative(const std::string &fname) const;

    /** Adds the mean value to the logged values */
    void log_mean_values(std::map<std::string, double> &logvalues) const;

    /** Get the lagrange multiplier */
    double get_lagrange_multiplier(unsigned int field, const MMSP::grid<dim, MMSP::vector<double> >&deriv_free) const;

    /** Add right hand side contribution from volume constraint */
    void add_volume_conservering_contribution(std::vector<double> &rhs, double lagrange, unsigned int field) const;
};

#endif