import unittest
import numpy as np
from apal import CahnHilliard
from apal_cxx import PyCahnHilliard


class TestCahnHilliard(unittest.TestCase):
    def test_fit(self):
        x = np.linspace(0.0, 1.0)
        a = 1.0
        b = 2.0
        c = 3.0
        d = 4.0
        e = 5.0
        G = a + b*x + c*x**2 + d*x**3 + e*x**4

        cahn = CahnHilliard(degree=4)
        cahn.fit(x, G)
        expected = [e, d, c, b, a]
        self.assertTrue(np.allclose(cahn.coeff, expected))

    def test_cython(self):
        coeff = [4.0, 3.0, 2.0, 1.0]
        cahn = CahnHilliard(degree=3, coeff=coeff)
        cython_cahn = PyCahnHilliard(coeff)

        x_values = [2.0, 3.0, -1.0]
        for x in x_values:
            self.assertAlmostEqual(cahn.evaluate(x), cython_cahn.evaluate(x))

    def test_raise_inconsistent_args(self):
        coeff = [4.0, 3.0, 2.0, 1.0]
        self.assertRaises(ValueError, CahnHilliard, coeff=coeff)

    def test_derivative(self):
        coeff = [4.0, 3.0, 2.0, 1.0]
        cahn = CahnHilliard(degree=3, coeff=coeff)
        cython_cahn = PyCahnHilliard(coeff)
        x_values = [2.0, 3.0, -1.0]
        for x in x_values:
            self.assertAlmostEqual(cahn.deriv(x), cython_cahn.deriv(x))

    def test_regularization_consistentcy(self):
        coeff = [4.0, 3.0, 2.0, 1.0]
        cahn = CahnHilliard(degree=3, coeff=coeff, bounds=[0.0, 1.0],
                            penalty=10.0, range_frac=0.1)
        cython_cahn = PyCahnHilliard(coeff, penalty=10.0, bounds=[0.0, 1.0],
                                     range_scale=0.1)

        x_values = [-0.1, 1.1, -0.5]
        for x in x_values:
            self.assertAlmostEqual(cahn.deriv(x), cython_cahn.deriv(x))

        for x in x_values:
            self.assertAlmostEqual(cahn.evaluate(x), cython_cahn.evaluate(x))

if __name__ == "__main__":
    unittest.main()
