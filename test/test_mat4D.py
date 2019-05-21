import unittest
import numpy as np
from apal_cxx import PyMat4D


class TestMat4D(unittest.TestCase):
    def test_numpy(self):
        array = np.random.rand(3, 3, 3, 3).astype(np.float64)
        mat4D = PyMat4D()
        mat4D.from_numpy(array)

if __name__ == "__main__":
    unittest.main()