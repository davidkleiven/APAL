from apal import TwoPhaseLandauPolynomialBase
from apal.landau_polynomial import array_func
import numpy as np

class SinglePrecipitatePoly(TwoPhaseLandauPolynomialBase):
    def __init__(self, c1=0.0, c2=1.0):
        TwoPhaseLandauPolynomialBase.__init__(self, c1=c1, c2=c2, conc_order2=2)
        self.conc_coeff1 = np.zeros(3)
        self.landau_amplitude = 1.0
    
    def _interpolation_function(self, x):
        return 3*x**2 - 2*x**3

    def _landau_barrier(self, x):
        return x**2 - 2*x**3 + x**4

    def _eval_phase1(self, conc):
        return np.polyval(self.conc_coeff1, conc)

    def _eval_phase2(self, conc):
        return np.polyval(self.conc_coeff2, conc)

    def _fit_second_order(self, x, y, x0=0.0):
        """
        Fit second order polynomial of the form
        A*(x - x0)^2 + B

        :param np.ndarray x: x-values
        :param np.ndarray y: y-values
        :param float x0: Position of the minimum
        """
        X = np.zeros((len(x), 2))
        X[:, 0] = 1.0
        X[:, 1] = (x-x0)**2
        coeff, _, _, _ = np.linalg.lstsq(X, y)
        return self._coeff2npy_poly(coeff, x0)

    def _coeff2npy_poly(self, coeff, x0):
        """
        Convert coefficients in the polynomial
        A*(x-x0)^2 + B to the form A*x^2 + B*x + C such that it
        can be used with np.polyval
        """
        out = np.zeros(3)
        out[0] = coeff[1]
        out[1] = -coeff[1]*x0
        out[2] = coeff[0]
        return out

    def _extract_subset(self, x, y, xlim=(0.0, 1.0)):
        """
        Extract subset of data with x-values given by xlim
        """
        indx_min = np.argmin(np.abs(x-xlim[0]))
        indx_max = np.argmin(np.abs(x-xlim[1]))

        if indx_max == len(y):
            indx_max = -1
        return x[indx_min:indx_max], y[indx_min:indx_max]

    @array_func
    def equil_shape_order(self, conc):
        dE = self._eval_phase2(conc) - self._eval_phase1(conc)
        W = self.landau_amplitude
        D = 1.5*((dE/W)**2 + dE/W + 1.5)

        if D < 0.0:
            return 0.0
        return 1.5*(dE/W + 1) + np.sqrt(D)

    def evaluate(self, conc, shape=None):
        if shape is None:
            shape = self.equil_shape_order(conc)
        return self._eval_phase1(conc)*(1 - self._interpolation_function(shape)) + \
            self._eval_phase2(conc)*self._interpolation_function(shape) + \
            self._landau_barrier(shape)*self.landau_amplitude

    def fit(self, x, G, lim1=(0.0, 0.1), lim2=(0.9, 1.0), show=False):
        """
        Fit concentration free energy to data
        """
        xfit, yfit = self._extract_subset(x, G, xlim=lim1)
        self.conc_coeff1 = self._fit_second_order(xfit, yfit, x0=self.c1)
        xfit, yfit = self._extract_subset(x, G, xlim=lim2)
        self.conc_coeff2 = self._fit_second_order(xfit, yfit, x0=self.c2)

        if show:
            from matplotlib import pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(x, G)
            ax.plot(x, self._eval_phase1(x))
            ax.plot(x, self._eval_phase2(x))

    def fit_landau_barrier(self, eta, G):
        # Shift the minimum of G to 0
        G += np.min(G)

        f_barrier = self._landau_barrier(eta)
        
        # Find the optimal barrier
        S1 = np.sum(f_barrier*G)
        S2 = np.sum(f_barrier**2)
        self.landau_amplitude = S1/S2
