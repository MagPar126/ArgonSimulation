"""
@authors: Magdalena and Johannes

Functions for calculating different forces and energies of the simulated system.

"""

import numpy as np

def LJ_force(diff_vectors):
    """
    Calculates the force of the Lennard-Jones potential on one particle from a given a number of relative positions.

    :param diff_vectors: List of difference vectors (np.array) for one particle (excluding the distance to the particle itself).
    :return: Force vector.
    """
    F = np.zeros(diff_vectors[0].shape)

    for dv in diff_vectors:
        r = np.linalg.norm(dv) #+ 1e-5
        F += dv/r * 4 * (12*r**(-13) - 6*r**(-7))

    return F

def kinetic_energy(velocity):
    """
    Calculates the classical kinetic energy of one particle from a given velocity.

    :param velocity: Velocity vector (numpy array).
    :return: Kinetic energy.
    """
    kin_energy = 0.5*(np.linalg.norm(velocity)**2)
    return kin_energy

def LJ_potential_energy(diff_vectors):
    """
    Calculates the Lennard-Jones potential energy of one particle from a given a number of relative positions.

    :param diff_vectors: List of difference vectors (np.array) for one particle (excluding the distance to the particle itself).
    :return: Potential energy.
    """
    pot_energy = 0
    for dv in diff_vectors:
        r = np.linalg.norm(dv)
        pot_energy += 4 * (r**(-12) - r**(-6))
    return pot_energy
    
