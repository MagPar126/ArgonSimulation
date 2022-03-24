"""
@authors: Magdalena and Johannes

Functions for calculating different forces and energies of the simulated system.

"""

import numpy as np

def LJ_force(diff_vectors):
    """
    Calculates the force of the Lennard-Jones potential on one particle from a given a number of relative positions.

    :param position: Current position of a particle.
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

def LJ_potential_energy(diff_vectors,position):
    """
    Calculates the Lennard-Jones potential energy of one particle from a given a number of relative positions.
    
    :param position: Current position of a particle.
    :param diff_vectors: List of difference vectors (np.array) for one particle (excluding the distance to the particle itself).
    :return: Potential energy.
    """
    pot_energy = 0
    for dv in diff_vectors:
        r = np.linalg.norm(dv)
        pot_energy += 4 * (r**(-12) - r**(-6))
    return pot_energy
    
def LJ_electric_field_force(electric_field,charge):
    """
    Takes a given electric field vector and charge assigned to the particle and created a force calculating function for LJ and electric potential.
    :param electric_field: Vector of constant electric field (list of floats).
    :param charge: Charge of a particle (float).
    :return: Function for calculating force, takes difference vectors, returns force vector.
    """
    def force_function(diff_vectors):
        F = LJ_force(diff_vectors)
        F+= np.array(electric_field, dtype=float) * charge
        return F
    return force_function

def LJ_electric_field_potential_energy(electric_field,charge):
    """
    Takes a given electric field vector and charge assigned to the particle and created a function calculating potential as a sum of LJ and electric potential.
    :param electric_field: Vector of constant electric field (list of floats).
    :param charge: Charge of a particle (float).
    :return: Function for potential, takes difference vectors and current position of the particle, returns potential energy.
    """
    def potential_function(diff_vectors,position):
        pot_energy = LJ_potential_energy(diff_vectors,position)
        for i in range(len(electric_field)):
            pot_energy -= position[i]*electric_field[i]

    
        return pot_energy
    return potential_function


    
    
