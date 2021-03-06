import unittest
import numpy as np
from itertools import product

from apal import Khachaturyan
from apal_cxx import PyKhachaturyan
from apal_cxx import pytest_functional_derivative
from apal_cxx import pytest_contract_tensors
from apal_cxx import pytest_B_tensor_element
from apal_cxx import pytest_strain_energy_sphere
from apal.tools import to_full_rank4


class TestKhacaturyan(unittest.TestCase):
    K = 50.0
    G = 26.0

    def get_isotropic_tensor(self):
        C = np.zeros((6, 6))
        C[0, 0] = C[1, 1] = C[2, 2] = self.K + 4*self.G/3.0
        C[0, 1] = C[0, 2] = \
        C[1, 0] = C[1, 2] = \
        C[2, 0] = C[2, 1] = self.K - 2.0*self.G/3.0
        C[3, 3] = C[4, 4] = C[5, 5] = 2*self.G
        return C

    @property
    def poisson(self):
        return 0.5*(3*self.K - 2*self.G)/(3*self.K + self.G)

    def isotropic_green_function(self, k):
        return np.eye(3)/self.G - 0.5*np.outer(k, k)/(self.G*(1.0 - self.poisson))

    def get_sphere_voxels(self, N):
        shape_func = np.zeros((N, N, N), dtype=np.uint8)
        indx = np.array(range(N))
        ix, iy, iz = np.meshgrid(indx, indx, indx)
        r_sq = (ix-N/2)**2 + (iy-N/2)**2 + (iz-N/2)**2
        r = N/8.0
        shape_func[r_sq<r] = 1
        return shape_func

    def eff_stress(self, elastic, misfit):
        return np.einsum('ijkl,kl->ij', elastic, misfit)

    def get_plate_voxels(self, N):
        shape_func = np.zeros((N, N, N), dtype=np.uint8)
        width = int(N/4)
        shape_func[:width, :width, :2] = 1
        return shape_func

    def get_needle_voxels(self, N):
        shape_func = np.zeros((N, N, N), dtype=np.uint8)
        width = int(N/4)
        shape_func[:width, :2, :2] = 1
        return shape_func

    def test_isotropic(self):
        misfit = np.eye(3)*0.05
        strain = Khachaturyan(elastic_tensor=self.get_isotropic_tensor(),
                              misfit_strain=misfit)
        k = np.array([5.0, -2.0, 7.0])
        khat = k/np.sqrt(k.dot(k))
        zeroth = strain.zeroth_order_green_function(khat)
        self.assertTrue(np.allclose(zeroth, self.isotropic_green_function(khat)))

    def eshelby_strain_energy_sphere(self, misfit):
        return 2*(1+self.poisson)*self.G*misfit**2/(1-self.poisson)

    def eshelby_strain_energy_plate(self, misfit):
        return 2*(1+self.poisson)*self.G*misfit**2/(1-self.poisson)

    def eshelby_strain_energy_needle(self, misfit):
        return 2*(1+self.poisson)*self.G*misfit**2/(1-self.poisson)

    def test_green_function_cpp(self):
        misfit = np.eye(3)*0.05
        elastic = to_full_rank4(self.get_isotropic_tensor())
        pykhach = PyKhachaturyan(3, elastic, misfit)
        k = np.array([-1.0, 3.0, 2.5])
        k /= np.sqrt(k.dot(k))
        gf = pykhach.green_function(k)
        self.assertTrue(np.allclose(gf, self.isotropic_green_function(k)))

    def test_frequency(self):
        misfit = np.eye(3)*0.05
        ft = np.zeros((8, 8, 8))
        elastic = to_full_rank4(self.get_isotropic_tensor())
        pykhach = PyKhachaturyan(3, elastic, misfit)

        freq = np.fft.fftfreq(ft.shape[0])
        for i in range(ft.shape[0]):
            indx = np.array([i, 0, 0])
            self.assertAlmostEqual(freq[i], pykhach.wave_vector(indx, ft.shape[0])[0])

    def test_sphere(self):
        eps = 0.05
        misfit = np.eye(3)*eps
        strain = Khachaturyan(elastic_tensor=self.get_isotropic_tensor(),
                              misfit_strain=misfit)
        sph = self.get_sphere_voxels(256)
        E = strain.strain_energy_voxels(sph)
        E_eshelby = self.eshelby_strain_energy_sphere(eps)
        self.assertAlmostEqual(E, E_eshelby, places=3)

    def test_sphere_pure_python(self):
        eps = 0.05
        misfit = np.eye(3)*eps
        strain = Khachaturyan(elastic_tensor=self.get_isotropic_tensor(),
                              misfit_strain=misfit)
        sph = self.get_sphere_voxels(32)
        E = strain.strain_energy_voxels(sph)
        E_eshelby = self.eshelby_strain_energy_sphere(eps)
        self.assertAlmostEqual(E, E_eshelby, places=3)

    def test_plate_voxels(self):
        eps = 0.05
        misfit = np.eye(3)*eps
        strain = Khachaturyan(elastic_tensor=self.get_isotropic_tensor(),
                              misfit_strain=misfit)
        plate = self.get_plate_voxels(256)
        E = strain.strain_energy_voxels(plate)
        E_eshelby = self.eshelby_strain_energy_plate(eps)
        self.assertAlmostEqual(E, E_eshelby, places=3)

    def test_needle_voxels(self):
        eps = 0.05
        misfit = np.eye(3)*eps
        strain = Khachaturyan(elastic_tensor=self.get_isotropic_tensor(),
                              misfit_strain=misfit)
        needle = self.get_needle_voxels(256)
        E = strain.strain_energy_voxels(needle)
        E_eshelby = self.eshelby_strain_energy_needle(eps)
        self.assertAlmostEqual(E, E_eshelby, places=3)

    def test_effective_stress(self):
        eps = 0.05
        misfit = np.eye(3)*eps
        misfit[0, 1] = 0.01
        misfit[0, 2] = -0.2
        misfit[1, 2] = 0.1
        misfit = 0.5*(misfit + misfit.T)

        elastic = to_full_rank4(self.get_isotropic_tensor())
        stress = self.eff_stress(elastic, misfit)

        khac = PyKhachaturyan(3, elastic, misfit)
        stress_cpp = khac.effective_stress()

        self.assertTrue(np.allclose(stress, stress_cpp))

    def test_contract_tensors(self):
        t1 = [[0.1, 0.2, 0.1],
              [0.1, 5.0, -0.2],
              [0.1, -0.2, -2.0]]
        t2 = [[-0.11, 2.0, 3.0],
              [2.0, 4.0, 0.2],
              [3.0, 0.2, -1.0]]
        cpp_contract = pytest_contract_tensors(t1, t2)

        pycontract = np.einsum("ij,ij", t1, t2)
        self.assertAlmostEqual(cpp_contract, pycontract)

    def test_B_tensor_element(self):
        gf = [[0.5, 0.2, 0.1],
              [0.2, -0.2, 0.3],
              [0.1, 0.3, 1.0]]
        t1 = [[0.1, 0.2, 0.1],
              [0.1, 5.0, -0.2],
              [0.1, -0.2, -2.0]]
        t2 = [[-0.11, 2.0, 3.0],
              [2.0, 4.0, 0.2],
              [3.0, 0.2, -1.0]]

        direction = np.array([0.3, 0.4, -0.1])
        direction /= np.sqrt(direction.dot(direction))

        cpp_element = pytest_B_tensor_element(direction, gf, t1, t2)

        py_elem = np.einsum("i,ij,jk,kl,l", direction, t1, gf, t2, direction)
        self.assertAlmostEqual(cpp_element, py_elem)

    def test_functional_derivative_one_field(self):
        elastic = to_full_rank4(self.get_isotropic_tensor())
        misfit = np.eye(3)
        misfit[0, 0] = 0.05
        misfit[1, 1] = -0.02
        misfit[2, 2] = 0.0

        init_field = np.zeros((128, 128))
        init_field[:15, :15] = 0.8

        try:
            result = pytest_functional_derivative(elastic, misfit, init_field.ravel())
        except RuntimeError as exc:
            # If fails, make sure that it is for the right reason
            self.assertTrue("The package was compiled without FFTW!" in str(exc))
            return

        func_deriv = result["func_deriv"]

        # Test 1 make sure that all entries outside 15x15 is zero
        self.assertTrue(np.allclose(init_field[15:, 15:], 0.0))

        init_field /= 0.8

        # Make sure field was passed correctly
        self.assertTrue(np.allclose(init_field, result["shape_squared_in"]))

        ft = np.fft.fft2(init_field**2)
        freq = np.fft.fftfreq(ft.shape[0])

        # Make sure that the real part of FFT match
        self.assertTrue(np.allclose(np.real(ft), result["ft_shape_real"]))

        V = 15*15
        stress = self.eff_stress(elastic, misfit)

        # Anlytical calculation
        for indx in product(range(ft.shape[0]), repeat=2):
            kvec = np.array([freq[indx[0]], freq[indx[1]], 0.0])
            k = np.sqrt(kvec.dot(kvec))

            if k < 1E-6:
                continue
            unit_vec = kvec/k
            G = self.isotropic_green_function(unit_vec)
            B = np.einsum('i,ij,jk,kl,l', unit_vec, stress, G, stress, unit_vec)
            ft[indx] *= B

        # Set the origin to the average value of the neighbours
        ft[0, 0] = 0.25*(ft[0, 1] + ft[1, 0] + ft[-1, 0] + ft[0, -1])
        self.assertTrue(np.allclose(np.real(ft), result["b_tensor_dot_ft_squared"]))

        ift_full = np.fft.ifft2(ft)

        self.assertTrue(np.allclose(np.imag(ift_full), 0.0))
        ift = np.real(ift_full)

        misfit_contrib = np.einsum('ijkl,ij,kl', elastic, misfit, misfit)*init_field**2
        self.assertTrue(np.allclose(misfit_contrib, result["misfit_energy_contrib"]))

        expect = 2*init_field*(misfit_contrib - ift)
        self.assertTrue(np.allclose(func_deriv, expect))

    def test_strain_field(self):
        misfit = np.zeros((3, 3))
        misfit[0, 0] = 0.05
        khach = Khachaturyan(elastic_tensor=self.get_isotropic_tensor(),
                             misfit_strain=misfit)
        
        shape = np.zeros((128, 128))
        shape[:, :20] = 1.0
        strain = khach.strain_field(shape)
        
        # Compare with exact solution
        self.assertTrue(np.allclose(strain[(0, 0)], 0.0))
        self.assertTrue(np.allclose(strain[(2, 2)], 0.0))
        self.assertTrue(np.allclose(strain[(0, 1)], 0.0))
        self.assertTrue(np.allclose(strain[(0, 2)], 0.0))
        self.assertTrue(np.allclose(strain[(1, 2)], 0.0))

        # TODO: Confirm that strain[(1, 1)] also satisfies
        # the solution. It seems to have the right structure
        # at least...

    def test_strain_energy(self):
        misfit = np.zeros((3, 3))
        misfit[0, 0] = misfit[1, 1] = misfit[2, 2] = 0.05
        elastic = to_full_rank4(self.get_isotropic_tensor())

        try:
            result = pytest_strain_energy_sphere(elastic, misfit)
        except RuntimeError as exc:
            # If fails, make sure that it is for the right reason
            self.assertTrue("The package was compiled without FFTW!" in str(exc))
            return
        

        # Expected misfit contribution
        vol = result["volume"]
        expected_misfit = 0.5*np.einsum("ijkl,ij,kl", elastic, misfit, misfit)
        self.assertAlmostEqual(expected_misfit, result["misfit_contrib"]/vol)
        

        eshelby_energy = self.eshelby_strain_energy_sphere(misfit[0, 0])
        self.assertAlmostEqual(eshelby_energy, result["energy"], places=2)

        
if __name__ == "__main__":
    unittest.main()

