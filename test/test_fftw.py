import unittest
from apal_cxx import pyfft1D, pyfft2D, pyfft3D
import numpy as np
import os

class TestFFTW(unittest.TestCase):
    def test_fft1D(self):
        os.environ['OMP_NUM_THREADS'] = "1"
        array = np.linspace(0.0, 10.0, 128)

        for direction in [-1, 1]:
            ft_real = pyfft1D(array.copy(), direction)

            if ft_real is None:
                self.skipTest("Could not perform 1D FFT")
            
            ft_npy = np.real(np.fft.fft(array))
            scale = ft_npy[0]/ft_real[0]
            self.assertTrue(np.allclose(ft_real*scale, ft_npy))

    def test_fft2D(self):
        os.environ['OMP_NUM_THREADS'] = "1"
        array = np.zeros((128, 128))
        array[40:70, 40:70] = 1.0

        # Use a symmetric array such that column or row major does not matter
        array += array.T
        array /= 2.0

        for direction in [-1, 1]:
            ft_real = pyfft2D(array.copy(), direction)

            if ft_real is None:
                self.skipTest("Coult not perform 2D FFT")
            
            ft_npy = np.real(np.fft.fft2(array))
            scale = ft_npy[0, 0]/ft_real[0, 0]
            ft_real *= scale
            self.assertTrue(np.allclose(ft_real, ft_npy))

    def test_fft3D(self):
        os.environ['OMP_NUM_THREADS'] = "1"
        array = np.zeros((128, 128, 128))
        array[40:70, 40:70, 40:70] = 1.0

        # Use a symmetric array such that column or row major does not matter
        array += array.T
        array /= 2.0

        for direction in [-1, 1]:
            ft_real = pyfft3D(array.copy(), direction)

            if ft_real is None:
                self.skipTest("Coult not perform 3D FFT")
            
            ft_npy = np.real(np.fft.fftn(array))
            scale = ft_npy[0, 0, 0]/ft_real[0, 0, 0]
            ft_real *= scale
            self.assertTrue(np.allclose(ft_real, ft_npy))

    def test_thread_consistency(self):
        import multiprocessing as mp
        num_cpu = mp.cpu_count()
        os.environ['OMP_NUM_THREADS'] = "1"
        array = np.linspace(0.0, 10.0, 128)
        ft_array = pyfft1D(array.copy(), 1)

        for num_threads in range(2, num_cpu+1):
            os.environ['OMP_NUM_THREADS'] = str(num_threads)
            ft_multithread = pyfft1D(array.copy(), 1)
            self.assertTrue(np.allclose(ft_array, ft_multithread))



if __name__ == "__main__":
    unittest.main()
