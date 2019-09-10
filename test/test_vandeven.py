import unittest
from apal_cxx import pytest_order1, pytest_order2


class TestVandeven(unittest.TestCase):
    def test_order1(self):
        self.assertTrue(pytest_order1())

    def test_order2(self):
        self.assertTrue(pytest_order2())

if __name__ == '__main__':
    unittest.main()