import unittest
from apal_cxx import pytest_omega_min
from apal_cxx import pytest_omega_max
from apal_cxx import pytest_eval_low_freq
from apal_cxx import pytest_eval_at_cut
from apal_cxx import pytest_eval_large_freq

class TestRaisedCosineFilter(unittest.TestCase):
    def test_omega_min(self):
        self.assertTrue(pytest_omega_min())
    
    def test_omega_max(self):
        self.assertTrue(pytest_omega_max())

    def test_eval_low_freq(self):
        self.assertTrue(pytest_eval_low_freq())

    def test_eval_at_cut(self):
        self.assertTrue(pytest_eval_at_cut())

    def test_eval_large_freq(self):
        self.assertTrue(pytest_eval_large_freq())


if __name__ == '__main__':
    unittest.main()