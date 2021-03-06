# distutils: language=c++

from apal.cython.regression_kernel cimport RegressionKernel

cdef extern from "regression_kernels.hpp":
    cdef cppclass QuadraticKernel(RegressionKernel):
        QuadraticKernel(double width)

        double evaluate(double x)

        double deriv(double x)

        object to_dict()

        void from_dict(object dict_repr) except+