import unittest
from apal_cxx import pytest_biharmonic_matrix
from apal_cxx import pytest_laplacian_matrix3D
from apal_cxx import pytest_add_same_element
from apal_cxx import pytest_add_mixed_element
from scipy.ndimage.filters import laplace
from scipy.linalg import toeplitz
import numpy as np

class TestBuildMatrices(unittest.TestCase):
    def get_3D_laplacian_matrix(self, N):
        col1 = np.zeros(N**3)
        col1[0] = -6
        col1[1] = 1
        col1[-1] = 1
        col1[N] = 1
        col1[-N] = 1
        col1[N**2] = 1
        col1[-N**2] = 1
        return toeplitz(col1)

    def test_biharmonic3D(self):
        """
        Check that the 3D biharmonic tensor can be constructed by an iterated
        laplacian (albeit not with nessecarily the same order of the nodes)
        """
        mat = pytest_biharmonic_matrix()

        # Subtract of 1 on the diagonal to get only the contribution
        # from the biharmonic operator
        mat -= np.identity(mat.shape[0])
        mat *= 0.5  # C++ code returns 2*biharmonic operator
        mat = mat.astype(np.int32)

        lapl3D = self.get_3D_laplacian_matrix(16)
        biharm_3D = lapl3D.dot(lapl3D).astype(np.int32)

        flat_template = biharm_3D.ravel()
        flat_cpp = mat.ravel()
        flat_template.sort()
        flat_cpp.sort()
        self.assertTrue(np.allclose(flat_cpp, flat_template))

    def test_laplacian_matrix3D(self):
        """
        Check that the concstructed 3D laplacian contains exactly the same 
        elements as the one constructed form a Toeplitz matrix
        """
        mat = pytest_laplacian_matrix3D().astype(np.int32)
        mat -= np.identity(mat.shape[0], dtype=np.int32)
        mat = -mat

        lapl3D = self.get_3D_laplacian_matrix(16).astype(np.int32)

        flat_lap3D = lapl3D.ravel()
        flat_lap3D.sort()

        flat_mat = mat.ravel()
        flat_mat.sort()

        max_diff = np.max(np.abs(flat_lap3D - flat_mat))
        msg = 'Max dev. {}'.format(max_diff)
        self.assertTrue(np.allclose(flat_lap3D, flat_mat), msg=msg)

    def test_add_same_element(self):
        res = pytest_add_same_element()
        self.assertTrue(np.allclose(res[0], res[1]))

    def test_add_mixed_element(self):
        res = pytest_add_mixed_element()
        self.assertTrue(np.allclose(res[0], res[1]))



if __name__ == "__main__":
    unittest.main()

