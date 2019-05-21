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

    def test_chglrealspace(self):
        chgl = self.get_chgl()
        chgl.build2D()

if __name__ == "__main__":
    unittest.main()