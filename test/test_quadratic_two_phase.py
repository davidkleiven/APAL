import unittest
from apal_cxx import PyQuadraticTwoPhasePoly, PyPolynomialTerm, PyPolynomial


class TestQuadraticTwoPhase(unittest.TestCase):

    def _get_example_poly1D(self):
        poly = PyPolynomial(1)
        powers = [2, 1]
        terms = [PyPolynomialTerm([powers[0]]), PyPolynomialTerm([powers[1]])]
        coeff = [2.0, -1.0]
        for c, t in zip(coeff, terms):
            poly.add_term(c, t)
        return poly, coeff, powers

    def _get_example_poly2D(self):
        poly = PyPolynomial(2)
        powers = [[2, 1], [1, 1], [1, 3]]
        coeff = [2.0, 1.0, 4.0]
        for c, p in zip(coeff, powers):
            poly.add_term(c, PyPolynomialTerm(p))
        return poly, coeff, powers

    def test_1D(self):
        poly = PyPolynomial(1)
        terms = [PyPolynomialTerm([2]), PyPolynomialTerm([1]), PyPolynomialTerm([0])]
        coeff = [2.0, 1.0, -1.0]

        for c, t in zip(coeff, terms):
            poly.add_term(c, t)
        
        quad = PyQuadraticTwoPhasePoly()
        quad.set_poly_phase1(poly)
        quad.set_poly_phase2(poly)

        x = 2.5
        expected = coeff[0]*x**2 + coeff[1]*x + coeff[2]

        # Divide result by two, since we add the same polynomial 
        # two places
        self.assertEqual(expected, 0.5*quad.evaluate_vec([x]))

    def test_error_raise(self):
        quad = PyQuadraticTwoPhasePoly()

        # Case 1: Raise because poly1 is not set
        with self.assertRaises(ValueError):
            quad.in_valid_state()

        poly = PyPolynomial(1)
        quad.set_poly_phase1(poly)

        # Case 2: Raise because poly2 is not set
        with self.assertRaises(ValueError):
            quad.in_valid_state()

        quad.set_poly_phase2(poly)
        quad.in_valid_state()

    def test_deriv_conc(self):
        quad = PyQuadraticTwoPhasePoly()

        poly1, coeff1, power1 = self._get_example_poly1D()
        poly2, coeff2, power2 = self._get_example_poly2D()
        
        quad.set_poly_phase1(poly1)
        quad.set_poly_phase2(poly2)

        x = [1.0, 2.0]
        deriv = quad.partial_deriv_conc_vec(x)

        expect1 = 0.0
        for c, p in zip(coeff1, power1):
            if p == 0:
                continue
            expect1 += p*c*x[0]**(p-1)
        
        expect2 = 0.0
        # Conc derivative is with respect to the first
        for c, p in zip(coeff2, power2):
            if p[0] == 0:
                continue
            expect2 += c*p[0]*x[0]**(p[0]-1)*x[1]**p[1]
        
        self.assertAlmostEqual(deriv, expect1 + expect2)

    def test_deriv_shape(self):
        quad = PyQuadraticTwoPhasePoly()

        poly2, coeff2, power2 = self._get_example_poly2D()
        quad.set_poly_phase2(poly2)

        x = [1.0, 4.0]

        # Derivative with respect to the first auxillary field
        deriv = quad.partial_deriv_shape_vec(x, 0)

        expected = 0.0
        for c, p in zip(coeff2, power2):
            if p[1] == 0:
                continue
            expected += c*(x[0]**p[0])*p[1]*x[1]**(p[1] - 1)
        self.assertAlmostEqual(deriv, expected)

        

if __name__ == "__main__":
    unittest.main()



        
