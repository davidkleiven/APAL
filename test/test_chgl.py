import unittest
import os
import numpy as np
from apal_cxx import PyCHGL
from apal_cxx import PyTwoPhaseLandau
from apal_cxx import PyPolynomial
from apal_cxx import PyKernelRegressor
from apal_cxx import PyGaussianKernel

class TestCHGL(unittest.TestCase):
    prefix = "chgl"
    L = 32

    def get_chgl(self):
        dim = 2
        num_gl_fields = dim
        M = 1.0
        alpha = 0.1
        dt = 0.001
        gl_damping = M
        grad_coeff = [[alpha, 0.5*alpha], [0.5*alpha, alpha]]
        return PyCHGL(dim, self.L, self.prefix, num_gl_fields, M, alpha, dt, gl_damping, 
                      grad_coeff)

    def test_run(self):
        chgl = self.get_chgl()
        chgl.print_polynomial()
        chgl.random_initialization([0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
        with self.assertRaises(RuntimeError):
            chgl.run(5, 1000)

        # We set a free energy form
        landau = PyTwoPhaseLandau()

        # This polynomial has the wrong size
        # make sure that an exception is raised
        poly = PyPolynomial(2)
        regressor = PyKernelRegressor(0.0, 1.0)
        kernel = PyGaussianKernel(2.0)

        # Case 1: Fail because no kernel is set
        with self.assertRaises(ValueError):
            chgl.set_free_energy(landau)

        # Case 2: Fail because no polynomial is set
        regressor.set_kernel(kernel)
        landau.set_kernel_regressor(regressor)
        with self.assertRaises(ValueError):
            chgl.set_free_energy(landau)

        poly = PyPolynomial(3)
        landau.set_polynomial(poly)
        chgl.set_free_energy(landau)

        try:
            chgl.run(5, 1000)
        except RuntimeError as exc:
            # The only way run should raise a runtime error at this stage is
            # if FFTW is not installed
            self.assertTrue("CHGL requires FFTW!" in str(exc))

    def test_npy_array(self):
        chgl = self.get_chgl()
        array = [np.random.rand(self.L, self.L),
                 np.random.rand(self.L, self.L),
                 np.random.rand(self.L, self.L)]
        chgl.from_npy_array(array)
        from_chgl = chgl.to_npy_array()
        for arr1, arr2 in zip(array, from_chgl):
            self.assertTrue(np.allclose(arr1, arr2))

    def test_exceptions(self):
        chgl = self.get_chgl()
        # Test wrong number of fields
        array = [np.random.rand(self.L, self.L),
                 np.random.rand(self.L, self.L)]
        with self.assertRaises(ValueError):
            chgl.from_npy_array(array)

        # Wrong dimension on one of the fields
        array = [np.random.rand(self.L, self.L),
                 np.random.rand(self.L, self.L),
                 np.random.rand(self.L, 2*self.L)]
        with self.assertRaises(ValueError):
            chgl.from_npy_array(array)

        # Wrong dimension on one field
        array = [np.random.rand(self.L, self.L),
                 np.random.rand(self.L, self.L),
                 np.random.rand(self.L, 2*self.L, 2)]
        with self.assertRaises(ValueError):
            chgl.from_npy_array(array)

    def tearDown(self):
        super(TestCHGL, self).tearDown()
        try:
            os.remove("chgl00000001000.grid")
            os.remove("chgl.grid")
        except OSError:
            pass

        try:
            os.remove("chgl_trackvalues.csv")
        except OSError:
            pass


if __name__ == "__main__":
    unittest.main()
    