import numpy as np

def LJ_force(diff_vectors, sigma=3.405e-10, epsilon=165.4e-23):
    """
    Calculate the force of the Lennard-Jones potential given a number of relative positions.

    :param diff_vectors:
    :param sigma:
    :param epsilon:
    :return: Force vector.
    """
    F = np.zeros(diff_vectors[0].shape)

    for dv in diff_vectors:
        r = np.linalg.norm(dv)
        F += dv/r * 4 * epsilon * (sigma**(12) * 12 * r**(-13) - sigma**6 * 6 * r**(-7))

    return F
