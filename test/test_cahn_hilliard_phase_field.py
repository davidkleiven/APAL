import unittest
import numpy as np
import os
from apal_cxx import PyCahnHilliardPhaseField
from apal_cxx import PyCahnHilliard


class TestCahnHilliardPhaseField(unittest.TestCase):
    def test_run_without_errors(self):
        coeff = [5.0, 4.0, 3.0, 2.0, 1.0]
        free = PyCahnHilliard(coeff)

        # Initialize a small 2D calculation
        L = 64
        M = 1.0
        dt = 0.01
        alpha = 1.0
        sim = PyCahnHilliardPhaseField(2, L, "cahnhill", free, M, dt, alpha)
        sim.run(100, 20)

    def test_run_without_errors3D(self):
        coeff = [5.0, 4.0, 3.0, 2.0, 1.0]
        free = PyCahnHilliard(coeff)

        # Initialize a small 2D calculation
        L = 16
        M = 1.0
        dt = 0.01
        alpha = 1.0
        sim = PyCahnHilliardPhaseField(3, L, "cahnhill", free, M, dt, alpha)
        sim.run(10, 2)

    def test_numpy_init(self):
        coeff = [5.0, 4.0, 3.0, 2.0, 1.0]
        free = PyCahnHilliard(coeff)

        L = 64
        M = 1.0
        dt = 0.01
        alpha = 1.0
        sim = PyCahnHilliardPhaseField(2, L, "cahnhill", free, M, dt, alpha)

        array = np.random.rand(L, L)
        sim.from_npy_array(array)
        from_sim = sim.to_npy_array()
        self.assertTrue(np.allclose(from_sim, array))

        array2 = np.random.rand(L)
        with self.assertRaises(ValueError):
            sim.from_npy_array(array2)

        array3 = np.random.rand(2*L, L)
        with self.assertRaises(ValueError):
            sim.from_npy_array(array3)

    def tearDown(self):
        super(TestCahnHilliardPhaseField, self).tearDown()
        try:
            os.remove("cahnhill0.vti")
            os.remove("cahnhill20.vti")
            os.remove("cahnhill40.vti")
            os.remove("cahnhill60.vti")
            os.remove("cahnhill80.vti")
            os.remove("cahnhill_adaptive_time.csv")
        except OSError:
            pass

        try:
            os.remove("cahnhill_trackvalues.csv")
        except OSError:
            pass

if __name__ == "__main__":
    unittest.main()
