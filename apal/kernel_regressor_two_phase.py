from apal import TwoPhaseLandauPolynomialBase
from apal.landau_polynomial import array_func
from apal.phasefield_util import fit_kernel, heaviside
from apal_cxx import PyKernelRegressor, PyGaussianKernel
import numpy as np

class TwoPhasePolynomialRegressor(TwoPhaseLandauPolynomialBase):
    def __init__(self, c1=0.0, c2=1.0, num_dir=3, init_guess=None,
                 conc_order1=2, conc_order2=2):
        TwoPhaseLandauPolynomialBase.__init__(self, c1=c1, c2=c2, num_dir=num_dir,
                                              init_guess=init_guess, conc_order1=conc_order1,
                                              conc_order2=conc_order2)
        self._kernel = None
        self._phase_one_regressor = None
        self.discontinuity_conc = 0.0
        self.discontinuity_jump = 0.0

    @property
    def phase_one_regressor(self):
        if self._phase_one_regressor is None:
            raise RuntimeError("It appears like no regressor was set. Call the fit function!")
        return self._phase_one_regressor

    @phase_one_regressor.setter
    def phase_one_regressor(self, new_obj):
        self._phase_one_regressor = new_obj

    @property
    def kernel(self):
        if self._kernel is None:
            raise RuntimeError("No kernel has been set!")
        return self._kernel

    @kernel.setter
    def kernel(self, new_obj):
        self._kernel = new_obj

    def _eval_phase1(self, conc):
        """Evaluate regressor in phase 1."""
        if self._phase_one_regressor is None:
            return 0.0
        return self._phase_one_regressor.evaluate(conc)

    def _deriv_phase1(self, conc):
        return self.phase_one_regressor.deriv(conc)

    @array_func
    def eval_at_equil(self, conc):
        """Evaluate the free energy at equillibrium order.

        :param float conc: Concentration
        """

        return TwoPhaseLandauPolynomialBase.eval_at_equil(self, conc) - \
            self.discontinuity_jump*heaviside(conc - self.discontinuity_conc)

    @array_func
    def evaluate(self, conc, shape=None):
        """
        Evaluate the free energy polynomial

        :param float conc: Concentration
        :param shape list: List with the shape order parameters.
            If None, the shape order parameters are set to their
            equillibrium
        """
        if shape is None:
            return self.eval_at_equil(conc)
        return TwoPhaseLandauPolynomialBase.evaluate(self, conc, shape=shape) - \
              self.discontinuity_jump*heaviside(conc - self.discontinuity_conc)

    def fit(self, conc, free_energy, weights={},
            init_shape_coeff=None, kernel_width=0.2, num_kernels=30,
            show=True, width=0.1, smear_width=0, shape="auto",
            lamb=None):
        """Fit the free energy functional.

        :param numpy.ndarray conc2. Concentrations in the second phase
        :param numpy.ndarray F2: Free energy in the second phase
        :param dict weights: Dictionary with penalty weights associated
            with the shape order parameter deviating from the desired
            behaviour. The following weights are available
            eq_phase1: A cost term equal to sum(n_eq**2) for concentrations
                less than self.c2 is added to the fitting cost function
            eq_phase2: A cost term equal to (n_eq(c_min) - 1.0)**2 where
                c_min is the concentration where F2 takes its minimum value.
                The effect of this is to penalize solutions that has a shape
                order paremter very different from one.
        """
        if init_shape_coeff is not None:
            if len(init_shape_coeff) != 2:
                raise ValueError("init_shape_coeff have to be None or "
                                 "of a list of length 2")

        if (width is None):
            peak = np.argmax(free_energy)
            conc_max = conc[peak]
            width = self.c2 - conc_max
        indx_min = np.argmin(np.abs(conc - self.c2 + width))
        indx_max = np.argmin(np.abs(conc - self.c2 - width))

        F_fit = free_energy[indx_min:indx_max]
        conc_fit = conc[indx_min:indx_max]
        minimum = np.argmin(np.abs(conc - self.c2))
        free_energy -= free_energy[minimum]

        if shape == "auto":
            shape = 27*(free_energy[indx_min] - free_energy[minimum])
            assert shape > 0.0

        # Update the coefficients
        A = self._get_slope_parameter(-shape, shape, conc[indx_min])
        self.conc_coeff2[0] = A
        self.conc_coeff2[1] = -A*self.c2
        self.coeff_shape[0] = -shape
        self.coeff_shape[2] = shape

        # Subtract of the reminder that needs to be fitted with 
        # Kernel regression
        reminder = free_energy - self.evaluate(conc)

        from matplotlib import pyplot as plt
        plt.plot(conc, reminder)
        plt.plot(conc, self.evaluate(conc))
        plt.plot(conc, free_energy)
        plt.show()

        # Smear the reminder in the critical area
        if smear_width > 0:
            smear_data = reminder[indx_min-smear_width:indx_min+smear_width]
            deriv1 = reminder[indx_min-smear_width] - \
                reminder[indx_min-smear_width-1]
            deriv2 = reminder[indx_min+smear_width+1] - \
                reminder[indx_min+smear_width]

            matrix = np.zeros((4, 4))
            rhs = np.zeros(4)
            imax = smear_width
            imin = -smear_width

            for i in range(4):
                # Continuity at beginning
                matrix[0, i] = imin**(3-i)

                # Continuity at end
                matrix[1, i] = imax**(3-i)

                # Smooth at beginning
                if i < 3:
                    matrix[2, i] = (3-i)*imin**(3-i-1)
                    matrix[3, i] = (3-i)*imax**(3-i-1)
            rhs[0] = reminder[indx_min-smear_width]
            rhs[1] = reminder[indx_min+smear_width]
            rhs[2] = deriv1
            rhs[3] = deriv2
            coeff = np.linalg.solve(matrix, rhs)
            x = np.linspace(imin, imax, 2*smear_width)
            reminder[indx_min-smear_width:indx_min+smear_width] = \
                sum(coeff[i]*x**(3-i) for i in range(4))

        self.kernel = PyGaussianKernel(kernel_width)

        # Remove discontinuity
        diff = np.abs(reminder - np.roll(reminder, 1))
        disc_indx = np.argmax(diff[1:])
        self.discontinuity_jump = reminder[disc_indx] - reminder[disc_indx+1]
        self.discontinuity_conc = conc[disc_indx]
        reminder[disc_indx+1:] += self.discontinuity_jump


        self.phase_one_regressor = fit_kernel(
            x=conc, y=reminder, num_kernels=num_kernels, kernel=self.kernel,
            lamb=lamb, extrapolate='linear', extrap_range=0.5)

        if show:
            from matplotlib import pyplot as plt
            fig = plt.figure()
            ax = fig.add_subplot(1, 2, 1)
            n_eq = np.array(self.equil_shape_order(conc))
            ax.plot(conc, free_energy, label="Free")
            ax.plot(conc, self.evaluate(conc), label="Evaluated")
            ax.plot(conc, self.phase_one_regressor.evaluate(conc), label="regressor")
            ax.axvline(conc[indx_min], ls="--", color="gray")
            ax.axvline(conc[indx_max], ls="--", color="gray")
            ax.legend()
            ax2 = fig.add_subplot(1, 2, 2)
            ax2.plot(conc, n_eq)

            plt.show()

    def to_dict(self):
        data = TwoPhaseLandauPolynomialBase.to_dict(self)


        data["discontinuity_jump"] = self.discontinuity_jump
        data["discontinuity_conc"] = self.discontinuity_conc

        # Store the kernel regressor
        try:
            data["kernel_regressor"] = self.phase_one_regressor.to_dict()
            data["kernel"] = self.kernel.to_dict()
        except RuntimeError:
            pass
        return data
        