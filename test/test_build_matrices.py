import unittest
from apal_cxx import pytest_biharmonic_matrix
from apal_cxx import pytest_laplacian_matrix3D
from apal_cxx import pytest_get_small_biharmonic
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

        self.assertTrue(np.allclose(flat_lap3D, flat_mat))


    def test_biharmonic_semi_positive_definite(self):
        mat = pytest_get_small_biharmonic()
        mat -= np.identity(mat.shape[0])

        eigvals = np.linalg.eigvalsh(mat)
        self.assertTrue(np.allclose(mat, mat.T))
        self.assertTrue(np.all(eigvals >= 0))

if __name__ == "__main__":
    unittest.main()

