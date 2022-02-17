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
        diff_vectors = []
        for i, particle in enumerate(self.Particles):
            other_particles = self.Particles.copy()
            other_particles.remove(particle)
            loc_diff_vectors=[] #diff vectors for this one particular particle
            for other_particle in other_particles:
                copies_distances = [] #list of distances of all copies from the particle
                copies_diff_vectors=[] #list of diff vectors of all copies from the particle
                
                basis = [-1, 0, 1]
                dims = []
                for j in range(self.dimension):
                    dims.append(len(basis))

                image_boxes = np.zeros(dims)
                it = np.nditer(image_boxes, flags=['multi_index'])
                for x in it:
                    displacement=(np.array(it.multi_index)-1)*self.length #displacements of copies in respect to the original one
                    mirrored_pos = other_particle.pos + displacement #possition of the mirrored particle
                    diff_vector = particle.pos - mirrored_pos #diff vector of the mirrored particle
                    mirrored_distance = np.linalg.norm(diff_vector) 
                    copies_distances.append(mirrored_distance)
                    copies_diff_vectors.append(diff_vector)
                
                index_min_copy = np.argmin(copies_distances)
                loc_diff_vectors.append(copies_diff_vectors[index_min_copy])
            
            diff_vectors.append(loc_diff_vectors)
        return diff_vectors

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



