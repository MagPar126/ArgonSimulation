import numpy as np

def LJ_force(diff_vectors):
    """
    Calculate the force of the Lennard-Jones potential given a number of relative positions.

    :param diff_vectors:
    :return: Force vector.
    """
    F = np.zeros(diff_vectors[0].shape)

    for dv in diff_vectors:
        r = np.linalg.norm(dv) + 1e-5
        F += dv/r * 4 * (12*r**(-13) - 6*r**(-7))

    return F

def kinetic_energy(velocity):
    kin_energy = 0.5*(np.linalg.norm(velocity)**2)
    return kin_energy

def LJ_potential_energy(diff_vectors):
    pot_energy = 0
    for dv in diff_vectors:
        r = np.linalg.norm(dv)
        pot_energy += 4 * (r**(-12) - r**(-6))
    return pot_energy
    
