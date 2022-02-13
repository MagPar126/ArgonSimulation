import numpy as np

import forces as forces

class Particle:
    """
        Particle
    """

    def __init__(self, pos0, vel0, mass=1):
        # list for saving all positions of particle
        self.trajectory = [pos0]

        # current position and velocity
        self.pos = pos0
        self.vel = vel0
        self.mass = mass


class MolecularDynamics:
    """
    Whole simulation of n particles in a box interacting via Lennard-Jones-potential.
    """

    def __init__(self, length, num_particles, time_step, dimension=2, mass=1, potential='LJ'):
        self.length = length
        self.num_particles = num_particles
        self.h = time_step
        self.dimension = dimension
        self.mass = mass

        if potential == 'LJ':
            self.force = forces.LJ_force
        else:
            print("Wrong wording or force not implemented yet!")
            exit(-1)

        self.initialization()

    def initialization(self):
        """
        Initialize N particles at random positions with random velocities.

        :param num_particles:
        :param mass:
        :return: 0 if success
        """

        self.Particles = []

        np.random.seed(0)
        positions = self.length * np.random.rand(self.num_particles, self.dimension)
        velocities = 2 * self.length/0.01 * np.random.rand(self.num_particles, self.dimension) - self.length/0.01 # <<<<<FIXME: include a meaningful upper boundary

        for i in range(self.num_particles):
            (self.Particles).append(Particle(positions[i,:], velocities[i,:], self.mass))
            # print("v: {}".format(self.Particles[i].pos))

        return 0

    def minimum_image_convention(self):
        """
        Minimum image convention for all particles.

        :return: List of list of difference vectors of the (mirrored) particles.
        """

        for i, particle in enumerate(self.Particles):
            other_particles = self.Particles.copy()
            other_particles.remove(particle)

            for other_particle in other_particles:
                mirrored_positions = [other_particle.pos]

                # num_diff = 3**self.dimension
                basis = self.length * [-1, 0, 1]

                #comb_basis = []
                #for j in range(self.dimension):
                #    comb_basis.append(basis)
#
                #basis_x, basis_y = np.meshgrid(basis, basis)
                #
                #new_diff = basis_x + basis_y


                #for pos in mirrored_positions:
                #    if len(mirrored_positions) < 3**self.dimension:
                #        for j in range(self.dimension):
                #            for de in basis:
                #                ev = np.zeros(self.dimension)
                #                ev[j] = 1
                #                new_pos = pos + de*ev
                #                print("pos {}".format(pos))
                #                print(de)
                #                print(ev)
                #                print("new pos: {}".format(new_pos))
                #                mirrored_positions.append(new_pos)
                #                print(mirrored_positions)
                #    else:
                #        break
#
#
                #print(mirrored_positions)

        return 0

    def simulation_step(self):
        list_diff_vectors = self.minimum_image_convention()

        for i, particle in enumerate(self.Particles):
            diff_vectors = list_diff_vectors[i]

            particle.pos += particle.vel * self.h # <<<<<FIXME: Save the position IN THE BOX, and not outside
            particle.pos = particle.pos % self.length # <<<<<<<<< checkme
            particle.vel += self.h * self.force(diff_vectors) / self.mass

            (particle.trajectory).append(particle.pos)

        return 0

    def simulate(self, num_steps):

        for t in range(num_steps):
            self.simulate_step()

        self.save_trajectories()

        return 0

    def save_trajectories(self):
        pass



