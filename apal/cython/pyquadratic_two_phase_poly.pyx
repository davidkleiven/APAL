from cython.operator cimport dereference

cdef extern from "polynomial.hpp":
    cdef cppclass Polynomial:
        pass

cdef extern from "quadratic_two_phase_poly.hpp":
    cdef cppclass QuadraticTwoPhase:
        QuadraticTwoPhase()

        double evaluate_vec(double x[])

        double partial_deriv_conc_vec(double x[])

        double partial_deriv_shape_vec(double x[], unsigned int dir)

        void in_valid_state() except+

        void set_poly_phase1(const Polynomial &poly) except+

        void set_poly_phase2(const Polynomial &poly)


cdef class PyQuadraticTwoPhasePoly:
    cdef QuadraticTwoPhase *thisptr

    def __cinit__(self):
        self.thisptr = new QuadraticTwoPhase()

    def __dealloc__(self):
        del self.thisptr

    def evaluate_vec(self, x):
        cdef double x_vec[10]

        for i in range(len(x)):
            x_vec[i] = x[i]
        return self.thisptr.evaluate_vec(x_vec)

    def partial_deriv_conc_vec(self, x):
        cdef double x_vec[10]
        
        for i in range(len(x)):
            x_vec[i] = x[i]

        return self.thisptr.partial_deriv_conc_vec(x_vec)

    def partial_deriv_shape_vec(self, x, direction):
        cdef double x_vec[10]

        for i in range(len(x)):
            x_vec[i] = x[i]
        return self.thisptr.partial_deriv_shape_vec(x_vec, direction)

    def in_valid_state(self):
        self.thisptr.in_valid_state()

    def set_poly_phase1(self, PyPolynomial poly):
        self.thisptr.set_poly_phase1(dereference(poly.thisptr))

    def set_poly_phase2(self, PyPolynomial poly):
        self.thisptr.set_poly_phase2(dereference(poly.thisptr))

