Simple Cahn Hilliard Model
===========================

This short tutorial shows how one can use APAL to run a simple Cahn-Hilliard calculation.
First we need to import a few Cython-wrapped classes

>>> from apal_cxx import PyCahnHilliard  # Nessecary to hold the information about the free energy
>>> from apal_cxx import PyCahnHilliardPhaseField  # Actually perform the simulation
>>> import numpy as np  # Numpy is used to initialize the field

Let's define all parameters nessecary to solve the Cahn-Hilliard equation

>>> prefix = "./simple_cahn_hilliard"  # Prefix inserted in front of all outputs
>>> L = 64       # Simulation domain size
>>> alpha = 0.5  # Gradient coefficient
>>> dim = 2      # Number of dimensions
>>> M = 0.1      # Mobility
>>> dt = 0.01    # Initial timestep

Now we need to define the functional form of the free energy. Here, we use a simple double well potential

:math:`f(c) = (1 - c^2)^2 = 1 - 2c^2 + c^4`

which is represented as an array of coefficients

>>> coeff = [1, 0, -2, 0, 1]
>>> free_energy = PyCahnHilliard(coeff)

Now we initialise the Cython class that is going to carry out the calculation

>>> cahn_hilliard = PyCahnHilliardPhaseField(2, L, prefix, free_energy, M, dt, alpha)
>>> min_stepsize = 1E-10  # Terminate if the timestep becomes lower than this
>>> max_change_per_step = 0.05  # Maximum allowed change per step (otherwise the timestep will be refined)
>>> cahn_hilliard.set_adaptive(min_stepsize, max_change_per_step)

Now we prepare build the nessecary matrices in order to propagate the equations using a semi-implicit scheme

>>> cahn_hilliard.build2D()  # If this was 3D, call build3D() instead

Then we create the initial concentration field

>>> init_conc = np.random.rand(L, L)
>>> cahn_hilliard.from_npy_array(init_conc)

Finally, we can solve the equation

>>> total_num_steps = 10  # The total number of timesteps
>>> output_every = 4      # Output a grid file every fourth timestep
>>> cahn_hilliard.run(total_num_steps, output_every) #doctest:+SKIP

