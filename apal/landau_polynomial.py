import numpy as np
from scipy.optimize import minimize, LinearConstraint, fsolve
from scipy.optimize import NonlinearConstraint
import scipy
import time

SCIPY_VERSION = scipy.__version__


def array_func(func):
    """
    Decorator to handle various array arguments
    """
    def unwrap(*args, **kwargs):
        self = args[0]
        conc = args[1]
        shape = kwargs.pop("shape", None)

        if isinstance(conc, list):
            conc = np.array(conc)
        if isinstance(shape, list):
            shape = np.array(shape)

        if np.isscalar(conc) and shape is None:
            return func(self, conc, **kwargs)

        elif isinstance(conc, np.ndarray) and shape is None:
            if len(conc.shape) != 1:
                raise ValueError("Concentration has to be a 1D array!")
            return [func(self, conc[i], **kwargs)
                    for i in range(conc.shape[0])]

        elif np.isscalar(conc) and np.isscalar(shape):
            return func(self, conc, shape=[shape, 0.0, 0.0], **kwargs)

        elif np.isscalar(conc) and isinstance(shape, np.ndarray):
            if len(shape.shape) == 1 and shape.shape[0] == 3:
                return func(self, conc, shape=shape, **kwargs)
            elif len(shape.shape) == 2 and shape.shape[1] == 3:
                return [func(self, conc, shape[i, :], **kwargs)
                        for i in range(shape.shape[0])]
            else:
                raise ValueError("When shape is a Numpy array it has to be "
                                 "either of length 3 or of length Nx3! "
                                 "Got: {}".format(shape.shape))

        elif isinstance(conc, np.ndarray) and isinstance(shape, np.ndarray):
            if conc.shape[0] != shape.shape[0]:
                raise ValueError("The number entries in the shape array has to"
                                 " match the number of entries in the conc "
                                 "array!")
            if len(shape.shape) == 1:
                return [func(self, conc[i], shape=[shape[i], 0.0, 0.0],
                             **kwargs) for i in range(conc.shape[0])]
            elif shape.shape[1] == 3:
                return [func(self, conc[i], shape=shape[i, :], **kwargs)
                        for i in range(conc.shape[0])]
            else:
                raise ValueError("Dimension of shape argument has to be either"
                                 "Nx3 or 3")
        else:
            raise ValueError("Concentation and shape arguments has to be "
                             "floats or arrays!")
    return unwrap


class TwoPhaseLandauPolynomialBase(object):
    """Class for fitting a Landau polynomial to free energy data

    In general terms a two phase landau polynomial is:
    
    1. A multidimensional object where there can be be up to
       three auxillary fields, but only one free variable
       At equillibrium, the auxillary fields are slaved by
       the concentration variable
    2. The "transition" where auxillary fields changed from one 
       value to another, is determined by a linear term of the
       form A*(x - x_c), where x_c is the transition point of
       the free variable.

    From the developers side, the idea behind having a base class is that
    it allows for different representations of the functional form of the 
    free variable, different fitting algorithms etc.

    :param float c1: Center concentration for the first phase
    :param float c2: Center concentration for the second phase
    :param np.ndarray init_guess: Initial guess for the parameters
        The polynomial fitting is of the form
        A*(x - c1)^2 + B*(x-c2)*y^2 + C*y^4 + D*y^6
        This array should therefore contain initial guess
        for the four parameter A, B, C and D.
    :param int conc_order1: Order of the polynomial in the first phase
    :param int conc_order2: Order of the polynomial in the second phase
    """
    def __init__(self, c1=0.0, c2=1.0, num_dir=3, init_guess=None,
                 conc_order1=2, conc_order2=2):
        self.conc_coeff2 = np.zeros(conc_order2+1)
        self.coeff_shape = np.zeros(5)
        self.conc_order1 = conc_order1
        self.conc_order2 = conc_order2
        self.c1 = c1
        self.c2 = c2
        self.init_guess = init_guess
        self.num_dir = num_dir
        self.boundary_coeff = None

    @array_func
    def equil_shape_order(self, conc):
        """Calculate the equillibrium shape concentration.

        The equillibrium shape order parameter is determined by finding the
        minima of the free energy curve at a given concentration. In case of
        multiple order parameters the value returned corresponds to a minima
        where all other shape order parameters are zero.

        :param float conc: Concentration
        """
        C = self.coeff_shape[0]
        D = self.coeff_shape[2]

        if abs(D) < 1E-8:
            n_eq = -0.5*self._eval_phase2(conc)/C
            if n_eq < 0.0:
                return 0.0
            return n_eq

        delta = (C/(3.0*D))**2 - \
            self._eval_phase2(conc)/(3.0*D)
        if delta < 0.0:
            return 0.0
        n_eq = -C/(3.0*D) + np.sqrt(delta)
        if n_eq < 0.0:
            return 0.0
        return np.sqrt(n_eq)

    @array_func
    def equil_shape_order_derivative(self, conc):
        """Calculate the partial derivative of the equillibrium
            shape parameter with respect to the concentration.

        NOTE: This return the derivative of the square of of the
            order parameter with respect to the concentration.
        """

        C = self.coeff_shape[0]
        D = self.coeff_shape[2]

        delta = (C/(3.0*D))**2 - \
            self._eval_phase2(conc)/(3.0*D)

        if delta < 0.0:
            return 0.0
        n_eq = self.equil_shape_order(conc)
        if n_eq <= 0.0:
            return 0.0
        return -0.5*self._deriv_phase2(conc) / \
            (3*np.sqrt(delta)*D*2*n_eq)

    def _eval_phase2(self, conc):
        """Evaluate the polynomial in phase2.

        :param float conc:
        """
        return np.polyval(self.conc_coeff2, conc)

    def _eval_phase1(self, conc):
        """Evaluate regressor in phase 1."""
        raise NotImplementedError("Has to be implemented in child classes")

    def _deriv_phase2(self, conc):
        """Evaluate the derivative in the second phase."""
        p2der = np.polyder(self.conc_coeff2)
        return np.polyval(p2der, conc)

    def _deriv_phase1(self, conc):
        raise NotImplementedError("Has to be implemented in child classes")

    @array_func
    def eval_at_equil(self, conc):
        """Evaluate the free energy at equillibrium order.

        :param float conc: Concentration
        """
        n_eq = self.equil_shape_order(conc)
        return self._eval_phase1(conc) + \
            self._eval_phase2(conc)*n_eq**2 + \
            self.coeff_shape[0]*n_eq**4 + \
            self.coeff_shape[2]*n_eq**6

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

        full_shape = np.zeros(3)
        full_shape[:len(shape)] = shape
        shape = full_shape
        return self._eval_phase1(conc) + \
            self._eval_phase2(conc)*np.sum(shape**2) + \
            self.coeff_shape[0]*np.sum(shape**4) + \
            self.coeff_shape[1]*(shape[0]**2 * shape[1]**2 +
                                 shape[0]**2 * shape[2]**2 +
                                 shape[1]**2 * shape[2]**2) + \
            self.coeff_shape[2]*np.sum(shape**6) + \
            self.coeff_shape[3]*(shape[0]**4 * (shape[1]**2 + shape[2]**2) +
                                 shape[1]**4 * (shape[0]**2 + shape[2]**2) +
                                 shape[2]**4 * (shape[0]**2 + shape[1]**2)) + \
            self.coeff_shape[4]*np.prod(shape**2)

    @array_func
    def partial_derivative(self, conc, shape=None, var="conc", direction=0):
        """Return the partial derivative with respect to variable."""
        allowed_var = ["conc", "shape"]
        if var not in allowed_var:
            raise ValueError("Variable has to be one of {}".format(allowed_var))

        if shape is None:
            shape = np.array([np.sqrt(self.equil_shape_order(conc))])

        if isinstance(shape, list):
            shape = np.array(shape)

        try:
            _ = shape[0]
        except (TypeError, IndexError):
            # Shape was a scalar, convert to array
            shape = np.array([shape])

        full_shape = np.zeros(3)
        full_shape[:len(shape)] = shape
        shape = full_shape

        if var == "conc":
            p2_der = np.polyder(self.conc_coeff2)
            return self._deriv_phase1(conc) + self._deriv_phase2(conc)*np.sum(shape**2)

        elif var == "shape":
            d = direction
            return 2*self._eval_phase2(conc)*shape[d] + \
                4*self.coeff_shape[0]*shape[d]**3 + \
                2*self.coeff_shape[1]*shape[d]*(shape[(d+1) % 3] + shape[(d+2) % 3]) + \
                6*self.coeff_shape[2]*shape[d]**5 + \
                4*self.coeff_shape[3]*shape[d]**3*(shape[(d+1) % 3]**2 + shape[(d+2) % 3]**2) + \
                2*self.coeff_shape[3]*shape[d]*(shape[(d+1) % 3]**4 + shape[(d+2) % 3]**4) + \
                2*self.coeff_shape[4]*shape[d]*shape[(d+1) % 3]**2 * shape[(d+2) % 3]**2
        else:
            raise ValueError("Unknown derivative type!")

    def _get_slope_parameter(self, C, D, transition_conc):
        return C**2/(3*D*(transition_conc - self.c2))

    def fit(self, *args, **kwargs):
        raise NotImplementedError("Has to be implemented in child classes")

    def to_dict(self):
        """Store the required arguments that can be used to
            construct poly terms for phase field calculations."""
        from itertools import permutations
        data = {}
        data["terms"] = []

        num_terms = len(self.conc_coeff2)
        for power, c in enumerate(self.conc_coeff2.tolist()):
            for active_shape in range(1, 4):
                entry = {
                    "coeff": c,
                    "powers": [num_terms-power-1, 0, 0, 0]
                }
                entry["powers"][active_shape] = 2
                data["terms"].append(entry)

        power_templates = [
            [4, 0, 0],
            [2, 2, 0],
            [6, 0, 0],
            [4, 2, 0],
            [2, 2, 2]
        ]

        for i, p_template in enumerate(power_templates):
            used_perms = set()
            for perm in permutations(p_template):
                if perm in used_perms:
                    continue
                entry = {
                    "coeff": self.coeff_shape[i],
                    "powers": [0] + list(perm)
                }
                used_perms.add(perm)
                data["terms"].append(entry)

        return data

    def save_poly_terms(self, fname="pypolyterm.json"):
        import json
        with open(fname, 'w') as outfile:
            json.dump(self.to_dict(), outfile, indent=2)
        print("Coefficient stored in {}".format(fname))

    def _equil_shape_fixed_conc_and_shape_intermediates(self, conc, shape,
                                                        min_type):
        """Return helper quantities for the equillibrium shape."""
        K = self._eval_phase2(conc)
        K += self.coeff_shape[1]*shape**2
        K += self.coeff_shape[3]*shape**4

        Q = self.coeff_shape[0] + self.coeff_shape[3]*shape**2

        if min_type == "mixed":
            Q += 0.5*self.coeff_shape[1]
            Q += 0.5*self.coeff_shape[4]*shape**2

        D = 3.0*self.coeff_shape[2]
        return K, Q, D

    @array_func
    def equil_shape_fixed_conc_and_shape(self, conc, shape=None,
                                         min_type="pure"):
        """Return the equillibrium shape parameter.

        :param float conc: Concentration
        :param float shape: Shape parameter
        :param str min_type: Type of minimum. If pure, the third
            shape parameter is set to zero. If mixed, the two
            free shape parameters are required to be the same.
        """
        allowed_types = ["pure", "mixed"]

        if min_type not in allowed_types:
            raise ValueError("min_type has to be one of {}"
                             "".format(allowed_types))

        if shape is None:
            raise ValueError("Shape has to be passed!")

        shape = shape[0]
        K, Q, D = self._equil_shape_fixed_conc_and_shape_intermediates(
            conc, shape, min_type)

        delta = (Q/D)**2 - K/D

        if delta < 0.0:
            return 0.0
        n_sq = -Q/D + np.sqrt(delta)

        if n_sq < 0.0:
            return 0.0
        return np.sqrt(n_sq)

    @array_func
    def equil_shape_fixed_conc_and_shape_deriv(self, conc, shape=None,
                                               min_type="pure"):
        """Differentiate with respect to the fixed shap parameter."""
        if shape is None:
            raise ValueError("Shape has to be passed!")
        shape = shape[0]
        K, Q, D = self._equil_shape_fixed_conc_and_shape_intermediates(
            conc, shape, min_type)

        dQ_dn = 2*self.coeff_shape[3]*shape

        if min_type == "mixed":
            dQ_dn += 2*self.coeff_shape[4]*shape*0.5

        dK_dn = 2*self.coeff_shape[1]*shape + \
            4*self.coeff_shape[3]*shape**3

        n_eq = self.equil_shape_fixed_conc_and_shape(
            conc, shape=shape, min_type=min_type)

        if n_eq <= 0.0:
            return 0.0

        delta = (Q/D)**2 - K/D
        deriv = - dQ_dn/D + 0.5*(2*Q*dQ_dn/D**2 - dK_dn/D)/np.sqrt(delta)
        return 0.5*deriv/n_eq

    def plot_individual_polys(self):
        from matplotlib import pyplot as plt
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)
        conc = np.linspace(0.0, 1.0, 100)
        ph1 = self._eval_phase1(conc)
        ax.plot(conc, ph1, label="Phase1")
        ax.plot(conc, np.polyval(self.conc_coeff2, conc), label="Phase2")
        ax.legend()
        return fig

    def fit_fixed_conc_varying_eta(self, conc, eta, free_energy, weights={},
                                   constraints=[]):
        """Perform fit at fixed composition, but varying eta.

        :param float conc: Fixed concentration
        :param array eta: Array with eta values
        :param array free_energy: Free energy densities
        :param dict weights: Cost function to tune the fitting
            Possible constraints:
            w_peak: Penalize deviation between the peak of the predicted
                energy and the free_energy array
            center_peak: Penalize solutions where the peak is 
                positioned far from the center.
        """

        w_peak = weights.get("peak", 0.0)
        w_peak_at_center = weights.get("center_peak", 0.0)
        w_mixed_peaks = weights.get("mixed_peaks", 0.0)

        def mse_function(x):
            self.coeff_shape[1] = x[0]
            self.coeff_shape[3:] = x[1:]
            pred = np.array(self.evaluate(conc, shape=eta))

            pred = np.array(pred)
            mse = np.mean((pred - free_energy)**2)
            peak_dev = np.max(pred) - np.max(free_energy)

            cost = mse

            # for cnst in constraints:
            #     cost += cnst(self)
            return cost

        # Last term has to be positive
        num_coeff = len(self.coeff_shape) - 2
        x0 = np.zeros(num_coeff)
        A = np.zeros((4, num_coeff))
        ub = np.zeros(4)
        lb = np.zeros(4)
        ub[0] = np.inf
        A[0, -1] = 1.0

        # Monotonical increasing cross terms
        A[1, 0] = 1.0
        A[1, 1] = 3.0
        lb[1] = 0.0
        ub[1] = np.inf

        # Ensure function is monotonically increasing along the line
        # eta1 = eta2
        A[2, 0] = 1.0
        lb[2] = np.abs(self.coeff_shape[0])
        ub[2] = np.inf

        # Ensure function is monotonically increasing along the line
        # eta1 = eta2
        A[3, 1] = 1.0
        lb[3] = -np.abs(self.coeff_shape[0])
        ub[3] = np.inf

        cnst = LinearConstraint(A, lb, ub)
        cb = MinimizationProgressCallback(10, constraints=constraints)
        res = minimize(mse_function, x0, method="SLSQP", constraints=cnst,
                       callback=cb)

        # Update the shape coefficients
        self.coeff_shape[1] = res["x"][0]
        self.coeff_shape[3:] = res["x"][1:]

    def _square_conc_scan(self, conc, nmin, nmax, num):
        """Perform a scan for many shape parameters."""
        from itertools import product
        F = np.zeros((num, num))
        n = np.linspace(nmin, nmax, num)
        for indx in product(range(num), repeat=2):
            F[indx] = self.evaluate(conc, shape=[n[indx[0]], n[indx[1]], 0.0])
        return F

    def _interior_weight(self, conc, nmax, num):
        scanned = self._square_conc_scan(conc, 0.0, nmax, num)

        # Locate minima along axis
        F = np.zeros(num)
        n = np.linspace(0.0, nmax, num)
        for i in range(num):
            F[i] = self.evaluate(conc, shape=[n[i], 0.0, 0.0])

        # We don't count the ones that are higher than the minimum
        # on the axis
        minval = np.min(F)
        scanned[scanned > minval] = minval
        return np.mean((scanned-minval)**2)

    def gradient_coefficient(self, alpha, gamma, conc, surf_form, gamma0=None):
        from scipy.optimize import newton
        deriv = np.array(self.equil_shape_order_derivative(conc))
        n_eq = np.array(self.equil_shape_order(conc))
        energy_at_zero = self.evaluate(conc, shape=np.zeros_like(conc))
        energy_at_n_eq = np.array([self.evaluate(c, shape=[n, 0, 0]) for c, n in zip(list(conc), list(n_eq))])
        n_eq_active = np.zeros(len(conc), dtype=np.uint8)
        n_eq_active[energy_at_n_eq<energy_at_zero] = 1
        n_eq_active[n_eq<0.1] = 0


        # The derivative is infinite at the transition. Remove unrealistic high contributions
        deriv[deriv>100*np.median(deriv)] = 0.0

        def eq(x):
            beta = x[0]**2  # Avoid negative values inside square root
            integrand = np.sqrt(surf_form)*np.sqrt(alpha + beta*deriv**2)
            integrand[n_eq_active==0] = 0.0
            interface_energy = 2*np.trapz(integrand, x=conc)
            return interface_energy - gamma

        def jac(x):
            beta = x[0]**2  # Avoid negative values inside square root
            integrand = 0.5*np.sqrt(surf_form)/np.sqrt(alpha + beta*deriv**2)
            integrand[n_eq_active==0] = 0.0
            integrand *= deriv**2
            fprime = 2*np.trapz(integrand, x=conc)
            array_fprime = np.zeros_like(x)

            # We need to differentiate with respect to x, not beta
            # Apply kernel rule
            array_fprime[0] = 2*x[0]*fprime
            return array_fprime

        if gamma0 is None:
            gamma0 = np.sqrt(alpha)
        root, infodict, ier, mesg = fsolve(eq, gamma0, fprime=jac, full_output=True)
        beta = np.sqrt(root)
        print("Target surface tension: {}. Difference: {}".format(gamma, infodict["fvec"]))
        return beta

    def conc_grad_param(self, gamma, x, surf_form):
        """Obtain gradient coefficients via a linearised approximation to the governing equations.

        The following system of equations is solved

        gamma1 = 2*integral sqrt(surf_form)*sqrt(b1 + b2*(d_eta/dc)**2)
        gamma2 = 2*integral sqrt(surf_form)*sqrt(b1 + f*b2*(d_eta/dc)**2)
        f = gamma2/gamma1

        by expanding the last square root.

        :param float gamma1: First surface tension
        :param float gamma2: Second surface tension
        :param x np.ndarray: 1D array with concentrations
        :param surf_form np.ndarray: 1D array with surface formation energy
        """

        sqrt_surf = np.sqrt(surf_form)
        sqrt_alpha = 0.5*gamma/np.trapz(sqrt_surf, x=x)
        alpha = sqrt_alpha**2
        return alpha


class MinimizationProgressCallback(object):
    def __init__(self, log_every_sec, constraints=[]):
        self.log_every_sec = log_every_sec
        self.now = time.time()
        self.constraints = constraints

    def log(self, msg):
        print(msg)

    def __call__(self, xk):
        if time.time() - self.now > self.log_every_sec:
            for cnst in self.constraints:
                self.log(cnst.status_msg())
            self.now = time.time()