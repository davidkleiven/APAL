# distutils: language = c++

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "two_phase_landau_base.hpp":
    cdef cppclass TwoPhaseLandauBase:
        pass

cdef extern from "khachaturyan.hpp":
    cdef cppclass Khachaturyan:
        pass

cdef extern from "chgl.hpp":
    cdef cppclass CHGL[T]:
        CHGL(int L, string &prefix, unsigned int num_gl_fields, \
             double M, double alpha, double dt, double gl_damping, 
             const vector[vector[double]] &interface) except +

        void update(int steps)

        void run(unsigned int start, unsigned int nsteps, int increment) except+

        void random_initialization(unsigned int field, double lower, double upper)

        void from_file(string fname)

        void from_npy_array(object fields) except +

        object to_npy_array() except +

        void set_free_energy(TwoPhaseLandauBase *polyter) except+

        void print_polynomial()

        void save_free_energy_map(string &fname)

        void use_HeLiuTang_stabilizer(double coeff)

        void use_adaptive_stepping(double min_dt, int increase_dt_every, double low_en_cut)

        void set_filter(double width)

        void set_raised_cosine_filter(double omega_cut, double roll_off)

        void set_gaussian_filter(double width)

        void add_strain_model(Khachaturyan model, int field)

        void conserve_volume(unsigned int field)
