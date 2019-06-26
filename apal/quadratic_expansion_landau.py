from apal import TwoPhaseLandauPolynomialBase
from apal.landau_polynomial import array_func
import numpy as np

class QuadraticExpansionLandau(TwoPhaseLandauPolynomialBase):
    def __init__(self, c1=0.0, c2=1.0, num_dir=3):
        TwoPhaseLandauPolynomialBase.__init__(self, c1=c1, c2=c2, num_dir=num_dir)
        self.coeff_phase1 = np.zeros(3)

    def _eval_phase1(self, conc):
        return np.polyval(self.coeff_phase1, conc)
    
    def _deriv_phase1(self, conc):
        pder = np.polyder(self.coeff_phase1)
        return np.polyval(pder, conc)

    @array_func
    def eval_at_equil(self, conc):
        evaluated = TwoPhaseLandauPolynomialBase.eval_at_equil(self, conc)

        # The above function that order parameter changes when it can
        # however, it can be that it is still energetically favorable
        # to not change. Thus, evaluate when forcing it to stay
        # zero
        eval_zero = TwoPhaseLandauPolynomialBase.evaluate(self, conc, shape=0.0)
        
        if eval_zero < evaluated:
            return eval_zero
        return evaluated

    def fit(self, conc, free_energy, end_phase1=0.1, show=False):
        """Fit energy polynomial.
        
        :param np.ndarray conc: Array with concentrations
        :param np.ndarray free_energy: Array with free energies
        :param float end_phase1: Quadratic polymial will be fitted including
            concentrations up to this value.
        """
        minimum = np.argmin(np.abs(conc - self.c2))
        #free_energy -= free_energy[minimum]

        x = conc[conc < end_phase1]
        F = free_energy[conc < end_phase1]
        
        X = np.zeros((len(x), 2))
        X[:, 0] = 1
        X[:, 1] = (x - self.c1)**2
        coeff, _, _, _ = np.linalg.lstsq(X, F)
        self.coeff_phase1[0] = coeff[1]
        self.coeff_phase1[1] = -2*coeff[1]*self.c1
        self.coeff_phase1[2] = coeff[0]

        #if shape == "auto":
        # indx = np.argmin(np.abs(conc - transition_conc))
        # shape = 27*free_energy[indx]
        # shape = 27*(self._eval_phase1(transition_conc) - free_energy[indx])/4
        # #shape = 27*self._eval_phase1(transition_conc)
        # assert shape > 0.0

        n_eq = np.sqrt(2.0/3.0) # If C = -D

        # Ensure local minimum at c2
        A = -self._deriv_phase1(self.c2)/n_eq**2
        self.conc_coeff2 = np.zeros(2)
        self.conc_coeff2[0] = A
        self.conc_coeff2[1] = -A*self.c2

        indx = np.argmin(np.abs(conc - self.c2))
        shape = (free_energy[indx] - self._eval_phase1(self.c2))/(n_eq**6 - n_eq**4)

        self.coeff_shape[0] = -shape
        self.coeff_shape[2] = shape

        if show:
            from matplotlib import pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(1, 2, 1)
            n_eq = np.array(self.equil_shape_order(conc))
            ax.plot(conc, free_energy, label="Free")
            ax.plot(conc, self.evaluate(conc), label="Evaluated")
            ax.plot(conc, self._eval_phase1(conc), label="Phase 1 poly")
            ax.legend()
            ax2 = fig.add_subplot(1, 2, 2)
            ax2.plot(conc, n_eq)

            plt.show()
        
    def to_dict(self):
        data = TwoPhaseLandauPolynomialBase.to_dict(self)
        data["conc_phase1"] = self.coeff_phase1.tolist()
        return data



