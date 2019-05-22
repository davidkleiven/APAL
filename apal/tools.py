import numpy as np
import json


def cahn_hilliard_surface_parameter(conc, free_energy, interfacial_energy):
    """Calculate the gradient parameter according to Cahn-Hilliard

        Free Energy of a Nonuniform System. I. Interfacial Free Energy
        Cahn, John W. "JW Cahn and JE Hilliard,
        J. Chem. Phys. 28, 258 (1958)."
        J. Chem. Phys. 28 (1958): 258.

    :param np.ndarray conc: Concentrations
    :param np.ndarray free_energy: Free energy difference
    :param float interfacial_energy: Interfacial energy in the same energy units
    """

    integral = np.trapz(np.sqrt(free_energy), x=conc)
    return (0.5*interfacial_energy/integral)**2


def get_polyterms(fname):
    """Parse JSON file and return list of PyPolyterms.

    :param str fname: JSON file with the parameters
    """
    from apal_cxx import PyPolynomialTerm

    with open(fname, 'r') as infile:
        data = json.load(infile)

    poly_terms = []
    coefficients = []
    for entry in data["terms"]:
        poly_terms.append(PyPolynomialTerm(entry["powers"]))
        coefficients.append(entry["coeff"])
    return coefficients, poly_terms


def surface_formation(conc, free_energy):
    """Extract the surface energy from a binary curve. Assume that there 
       is one local minima in the first half of the free_energy array and one 
       in the second half.

       :param conc np.ndarray: 1D array with concentrations
       :param free_energy np.ndarray: 1D array with free_energies
    """

    N = len(free_energy)
    first_half = free_energy[:int(N/2)]
    last_half = free_energy[int(N/2):]

    min_first = np.argmin(first_half)
    min_second = np.argmin(last_half)

    y1 = free_energy[min_first]
    y2 = free_energy[int(N/2) + min_second]
    x1 = conc[min_first]
    x2 = conc[int(N/2) + min_second]

    tangent_slope = (y2 - y1)/(x2 - x1)
    
    tangent = tangent_slope*(conc - x1) + y1
    surf_form = free_energy - tangent

    i1 = min_first
    i2 = int(N/2) + min_second
    return conc[i1:i2], surf_form[i1:i2]

