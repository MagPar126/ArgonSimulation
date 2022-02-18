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
