import unittest
import os
import numpy as np
from apal_cxx import PyCHGLRealSpace
from apal_cxx import PyTwoPhaseLandau
from apal_cxx import PyPolynomial
from apal_cxx import PyKernelRegressor
from apal_cxx import PyGaussianKernel


class TestCHGLRealSpace(unittest.TestCase):
    prefix = "chglrealspace"
    L = 32

    def get_chgl(self):
        dim = 2
        num_gl_fields = dim
        M = 1.0
        alpha = 0.1
        dt = 0.001
        gl_damping = M
        grad_coeff = [[alpha, 0.5*alpha], [0.5*alpha, alpha]]
        return PyCHGLRealSpace(dim, self.L, self.prefix, num_gl_fields, M,
                               alpha, dt, gl_damping, grad_coeff)

    def get_chgl3D(self):
        dim = 3
        num_gl_fields = dim
        M = 1.0
        alpha = 0.1
        dt = 0.001
        gl_damping = M
        grad_coeff = [[alpha, 0.5*alpha, 0.7*alpha], [0.5*alpha, 0.7*alpha, alpha], [0.7*alpha, alpha, 0.5*alpha]]
        return PyCHGLRealSpace(dim, 16, self.prefix, num_gl_fields, M,
                               alpha, dt, gl_damping, grad_coeff) 

    def test_chglrealspace(self):
        chgl = self.get_chgl()
        chgl.build2D()
        chgl.random_initialization([0.0, 0.0, 0.0], [1.0, 0.85, 0.85])

        # Initialize a regression kernel and a regressor
        kernel = PyGaussianKernel(2.0)
        regressor = PyKernelRegressor(0.0, 1.0)

        # Initialize a two phase landau polynomial
        landau = PyTwoPhaseLandau()

        # Transfer the kernel and the regressor
        regressor.set_kernel(kernel)
        landau.set_kernel_regressor(regressor)

        # Initialize 3D polynomial (concentration, shape1, shape2)
        poly = PyPolynomial(3)
        landau.set_polynomial(poly)
        chgl.set_free_energy(landau)

        chgl.run(5, 1000)

    def test_chglrealspace3D(self):
        chgl = self.get_chgl3D()
        chgl.build3D()

         # Initialize a regression kernel and a regressor
        kernel = PyGaussianKernel(2.0)
        regressor = PyKernelRegressor(0.0, 1.0)

        # Initialize a two phase landau polynomial
        landau = PyTwoPhaseLandau()

        # Transfer the kernel and the regressor
        regressor.set_kernel(kernel)
        landau.set_kernel_regressor(regressor)

        # Initialize 4D polynomial (concentration, shape1, shape2, shape3)
        poly = PyPolynomial(4)
        landau.set_polynomial(poly)
        chgl.set_free_energy(landau)

        chgl.run(5, 1000)

    def tearDown(self):
        super(TestCHGLRealSpace, self).tearDown()
        try:
            os.remove('chglrealspacetrack_values.csv')
        except Exception:
            pass
if __name__ == "__main__":
    unittest.main()