#ifndef CHGL_H
#define CHGL_H
#include "phase_field_simulation.hpp"
#include "two_phase_landau_base.hpp"
#include "fftw_mmsp.hpp"
#include "concentration_type.hpp"
#include <set>

#ifdef HAS_FFTW
    #include <complex>
    #include <fftw3.h>
#endif
#include "fftw_complex_placeholder.hpp"

#include "MMSP.grid.h"
#include "MMSP.vector.h"
#include "thermal_noise_generator.hpp"
#include "multidirectional_khachaturyan.hpp"
#include <vector>
#include <string>

typedef std::vector<std::vector<double> > interface_vec_t;

template<int dim>
using ft_grid_t = MMSP::grid<dim, MMSP::vector<fftw_complex> >;

class FourierDomainFilter;

template<int dim>
class CHGL: public PhaseFieldSimulation<dim>{
public:
    CHGL(int L, const std::string &prefix, unsigned int num_gl_fields, \
         double M, double alpha, double dt, double gl_damping, 
         const interface_vec_t &interface);

    virtual ~CHGL();

    /** Add a new free energy term to the model */
    void set_free_energy(const TwoPhaseLandauBase *poly);

    /** Return an array of the free energy */
    void save_free_energy_map(const std::string &fname) const;

    /** Print the polynomial used to screen */
    void print_polynomial() const;

    /** Use the stablization scheme by He Lio and Tang */
    void use_HeLiuTang_stabilizer(double coeff);

    /** Refine the time step if the energy increases */
    void use_adaptive_stepping(double min_dt, unsigned int inc_every, double lower_energy_change_cut);

    /** Use Gaussian filter to supress high modes */
    void set_filter(double width);

    /** Add thermal noise to the equations */
    void set_cook_noise(double amplitude);

    /** Save a noise replica */
    void save_noise_realization(const std::string &fname, unsigned int field) const;

    /** Update the timestep */
    void set_timestep(double new_dt);

    /** Add a strain model */
    void add_strain_model(const Khachaturyan &model, unsigned int field){khachaturyan.add_model(model, field);};

    /** Implement the update function */
    virtual void update(int nsteps) override;

    /** Calculate the energy of the system */
    virtual double energy() const;

    /** Set a filter to remove high frequency components */
    void set_filter(const FourierDomainFilter &filter){ft_filter = &filter;};

    /** Set a raised cosine filter */
    void set_raised_cosine_filter(double omega_cut, double roll_off);

    /** Set a gaussian filter */
    void set_gaussian_filter(double width);

    /** Set vandeven filter */
    void set_vandeven_filter(unsigned int order);

    /** Conserve volume of the order parameter squared */
    void conserve_volume(unsigned int gl_field);

    /** Set the concentration equation to Allen-Cahn*/
    void set_conc_type_allen_cahn();

    /** Set the concentration equation to Cahn-Hilliard*/
    void set_conc_type_cahn_hilliard();
protected:
    double M;
    double alpha;
    double dt;
    double gl_damping;
    double stab_coeff{0.0};
    bool adaptive_dt{false};
    double minimum_dt{1E-8};
    double old_energy{0.0};
    double filter_width{1.0};
    bool use_filter{false};
    bool old_energy_initialized{false};
    unsigned int increase_dt{100000};
    double lower_energy_cut{0.0};
    unsigned int update_counter{0};
    const FourierDomainFilter* ft_filter{nullptr};
    bool own_ft_filter_ptr{false};
    std::set<unsigned int> conserved_gl_fields;
    ConcentrationType_t conc_type{ConcentrationType_t::CAHN_HILLIARD};

    interface_vec_t interface;
    const TwoPhaseLandauBase *free_energy{nullptr};
    MMSP::grid<dim, MMSP::vector<fftw_complex> > *cmplx_grid_ptr{nullptr};

    MultidirectionalKhachaturyan khachaturyan;

    FFTW *fft{nullptr};
    std::vector<ThermalNoiseGenerator*> cook_noise;

    /** Check that the provided interfaces vector matches requirements */
    void check_interface_vector() const;
    void from_parent_grid();
    void to_parent_grid() const;
    bool is_initialized{false};

    /** Gaussian filter for high modes */
    double gaussian_filter_weight(double k) const;
    void calculate_strain_contrib(const ft_grid_t<dim> &grid_in, ft_grid_t<dim> &out);

    /** Derivative of the volume interpolating function */
    double deriv_vol_interpolating_function(double eta) const;

    /** Calculate the lagrange multiplier */
    double lagrange_multiplier(double rhs, double zeroth_vol_interp);

    /** Calculate the total variation */
    void total_variation(const ft_grid_t<dim> &grid_in, std::vector<double> &tv) const;

    /** Decreasing TV */
    bool tv_is_decreasing(const std::vector<double> &old_tv, const std::vector<double> &new_tv) const;
};
#endif