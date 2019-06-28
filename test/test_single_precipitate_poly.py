import unittest
import numpy as np
from apal import SinglePrecipitatePoly

class TestSinglePrecipitatePoly(unittest.TestCase):
    def test_fit(self):
        c1 = 0.05
        c2 = 0.95
        A1 = 1.3
        A2 = 1.5
        C1 = -0.1
        C2 = 0.4

        poly = SinglePrecipitatePoly(c1=c1, c2=c2)
        x = np.linspace(0.0, 1.0, 100)
        y = np.zeros(100)

        # Insert first
        y[:50] = A1*(x[:50] - c1)**2 + C1
        y[50:] = A2*(x[50:] - c2)**2 + C2

        poly.fit(x, y, lim1=(0.0, 0.3), lim2=(0.7, 1.0))

        expected1 = [A1, -A1*c1, C1]
        self.assertTrue(np.allclose(poly.conc_coeff1, expected1))

        expected2 = [A2, -A2*c2, C2]
        self.assertTrue(np.allclose(poly.conc_coeff2, expected2))


if __name__ == "__main__":
    unittest.main()