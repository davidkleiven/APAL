# distutils: language = c++

from apal.cython.chgl cimport CHGL
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

# Hack to support integer template arguments with cython
#cdef extern from *:
#    ctypedef int intParameter1 "1"
#    ctypedef int intParameter2 "2"
#    ctypedef int intParameter3 "3"

cdef extern from "two_phase_landau_base.hpp":
    cdef cppclass TwoPhaseLandauBase:
        pass

cdef class PyCHGL:
    cdef CHGL[intParameter1] *thisptr1D
    cdef CHGL[intParameter2] *thisptr2D
    cdef CHGL[intParameter3] *thisptr3D
    cdef int dim

    def __cinit__(self, dim, L, prefix, num_gl_fields, M, alpha, dt, gl_damping, gradient_coeff):
        cdef vector[vector[double]] interface
        cdef vector[double] inner_vec
        for i in range(dim):
            inner_vec.push_back(0.0)

        for item in gradient_coeff:
            for i in range(len(item)):
                inner_vec[i] = item[i]
            interface.push_back(inner_vec)

        if dim == 1:
            self.thisptr1D = new CHGL[intParameter1](L, prefix, num_gl_fields, M, alpha, dt, gl_damping, interface)
        elif dim == 2:
            self.thisptr2D = new CHGL[intParameter2](L, prefix, num_gl_fields, M, alpha, dt, gl_damping, interface)
        elif dim == 3:
            self.thisptr3D = new CHGL[intParameter3](L, prefix, num_gl_fields, M, alpha, dt, gl_damping, interface)

    def __init__(self, dim, L, prefix, num_gl_fields, M, alpha, dt, gl_damping, gradient_coeff):
        self.dim = dim

    def __dealloc__(self):
        del self.thisptr1D
        del self.thisptr2D
        del self.thisptr3D
        self.thisptr1D = NULL
        self.thisptr2D = NULL
        self.thisptr3D = NULL


    def run(self, nsteps, increment, start=0):
        if self.dim == 1:
            self.thisptr1D.run(start, nsteps, increment)
        elif self.dim == 2:
            self.thisptr2D.run(start, nsteps, increment)
        elif self.dim == 3:
            self.thisptr3D.run(start, nsteps, increment)

    def from_file(self, fname):
        if self.dim == 1:
            self.thisptr1D.from_file(fname)
        elif self.dim == 2:
            self.thisptr2D.from_file(fname)
        elif self.dim == 3:
            self.thisptr3D.from_file(fname)
        else:
            raise ValueError("Dimension has to be 1, 2, or 3")

    def random_initialization(self, lower, upper):
        for i in range(len(lower)):
            if self.dim == 1:
                self.thisptr1D.random_initialization(i, lower[i], upper[i])
            elif self.dim == 2:
                self.thisptr2D.random_initialization(i, lower[i], upper[i])
            elif self.dim == 3:
                self.thisptr3D.random_initialization(i, lower[i], upper[i])
            else:
                raise ValueError("Dimension has to be 1, 2, or 3")

    def from_npy_array(self, arrays):
        if self.dim == 1:
            self.thisptr1D.from_npy_array(arrays)
        elif self.dim == 2:
            self.thisptr2D.from_npy_array(arrays)
        elif self.dim == 3:
            self.thisptr3D.from_npy_array(arrays)

    def to_npy_array(self):
        if self.dim == 1:   
            return self.thisptr1D.to_npy_array()
        elif self.dim == 2:
            return self.thisptr2D.to_npy_array()
        elif self.dim == 3:
            return self.thisptr3D.to_npy_array()

    def set_free_energy(self, PyTwoPhaseLandau term):
        if self.dim == 1:
            self.thisptr1D.set_free_energy(<TwoPhaseLandauBase*> term.thisptr)
        elif self.dim == 2:
            self.thisptr2D.set_free_energy(<TwoPhaseLandauBase*> term.thisptr)
        elif self.dim == 3:
            self.thisptr3D.set_free_energy(<TwoPhaseLandauBase*> term.thisptr)

    def print_polynomial(self):
        if self.dim == 1:
            self.thisptr1D.print_polynomial()
        elif self.dim == 2:
            self.thisptr2D.print_polynomial()
        elif self.dim == 3:
            self.thisptr3D.print_polynomial()

    def save_free_energy_map(self, fname):
        if self.dim == 1:
            self.thisptr1D.save_free_energy_map(fname)
        elif self.dim == 2:
            self.thisptr2D.save_free_energy_map(fname)
        elif self.dim == 3:
            self.thisptr3D.save_free_energy_map(fname)

    def use_HeLiuTang_stabilizer(self, coeff):
        if self.dim == 1:
            self.thisptr1D.use_HeLiuTang_stabilizer(coeff)
        elif self.dim == 2:
            self.thisptr2D.use_HeLiuTang_stabilizer(coeff)
        elif self.dim == 3:
            self.thisptr3D.use_HeLiuTang_stabilizer(coeff)

    def use_adaptive_stepping(self, min_dt, increase_every_update, low_en_cut=0.0):
        if self.dim == 1:
            self.thisptr1D.use_adaptive_stepping(min_dt, increase_every_update, low_en_cut)
        elif self.dim == 2:
            self.thisptr2D.use_adaptive_stepping(min_dt, increase_every_update, low_en_cut)
        elif self.dim == 3:
            self.thisptr3D.use_adaptive_stepping(min_dt, increase_every_update, low_en_cut)

    def set_filter(self, width):
        if self.dim == 1:
            self.thisptr1D.set_filter(width)
        elif self.dim == 2:
            self.thisptr2D.set_filter(width)
        elif self.dim == 3:
            self.thisptr3D.set_filter(width)
            
    def set_free_energy_quadratic(self, PyQuadraticTwoPhasePoly term):
        if self.dim == 1:
            self.thisptr1D.set_free_energy(<TwoPhaseLandauBase*> term.thisptr)
        elif self.dim == 2:
            self.thisptr2D.set_free_energy(<TwoPhaseLandauBase*> term.thisptr)
        elif self.dim == 3:
            self.thisptr3D.set_free_energy(<TwoPhaseLandauBase*> term.thisptr)

    def set_raised_cosine_filter(self, omega_cut, roll_off):
        if self.dim == 1:
            self.thisptr1D.set_raised_cosine_filter(omega_cut, roll_off)
        elif self.dim == 2:
            self.thisptr2D.set_raised_cosine_filter(omega_cut, roll_off)
        elif self.dim == 3:
            self.thisptr3D.set_raised_cosine_filter(omega_cut, roll_off)

    def set_gaussian_filter(self, omega_cut):
        if self.dim == 1:
            self.thisptr1D.set_gaussian_filter(omega_cut)
        elif self.dim == 2:
            self.thisptr2D.set_gaussian_filter(omega_cut)
        elif self.dim == 3:
            self.thisptr3D.set_gaussian_filter(omega_cut)

    def set_vandeven_filter(self, order):
        if self.dim == 1:
            self.thisptr1D.set_vandeven_filter(order)
        elif self.dim == 2:
            self.thisptr2D.set_vandeven_filter(order)
        elif self.dim == 3:
            self.thisptr3D.set_vandeven_filter(order)

    def add_strain_model(self, PyKhachaturyan obj, field):
        if self.dim == 1:
            self.thisptr1D.add_strain_model(deref(obj.thisptr), field)
        elif self.dim == 2:
            self.thisptr2D.add_strain_model(deref(obj.thisptr), field)
        elif self.dim == 3:
            self.thisptr3D.add_strain_model(deref(obj.thisptr), field)

    def conserve_volume(self, gl_field):
        if self.dim == 1:
            self.thisptr1D.conserve_volume(gl_field)
        elif self.dim == 2:
            self.thisptr2D.conserve_volume(gl_field)
        elif self.dim == 3:
            self.thisptr3D.conserve_volume(gl_field)

    def set_conc_type_cahn_hilliard(self):
        if self.dim == 1:
            self.thisptr1D.set_conc_type_cahn_hilliard()
        elif self.dim == 2:
            self.thisptr2D.set_conc_type_cahn_hilliard()
        elif self.dim == 3:
            self.thisptr3D.set_conc_type_cahn_hilliard()

    def set_conc_type_allen_cahn(self):
        if self.dim == 1:
            self.thisptr1D.set_conc_type_allen_cahn()
        elif self.dim == 2:
            self.thisptr2D.set_conc_type_allen_cahn()
        elif self.dim == 3:
            self.thisptr3D.set_conc_type_allen_cahn()