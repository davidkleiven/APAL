import unittest
from apal_cxx import pytest_biharmonic_matrix
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
        

        
       

        # # Create an 3D image
        # image = np.zeros((16, 16, 16))

        # L = mat.shape[0]
        # size = int(L/4)

        # # Insert sime features
        # image[:size,:size:size] = 1.0
        # image[2*size:, :, :] = 1.0

        # # Flattened image
        # flattened = image.copy()
        # flattened = flattened.ravel()

        # # Biharmonic
        # biharm_filter = 0.5*mat.dot(flattened)
        # biharm_filter = np.reshape(biharm_filter, image.shape)

        # # First make sure that the filtered image is different
        # # from the original
        # self.assertFalse(np.allclose(image, biharm_filter))

        # # Apply laplacian filter
        # output = image.copy()
        # laplace(image, output, mode='wrap')
        # image = output.copy()

        # # Again make sure that it still is different
        # self.assertFalse(np.allclose(image, biharm_filter))

        # # Apply laplacian filter again
        # laplace(image, output, mode='wrap')
        # image = output.copy()

        # # Now it should be the same has the biharmonic filtered
        # #print(image)
        # print(np.max(np.abs(biharm_filter)))
        
        # self.assertTrue(np.allclose(image, biharm_filter))


if __name__ == "__main__":
    unittest.main()

