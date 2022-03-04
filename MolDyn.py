import numpy as np

import forces as forces

import plotting as plotting

import matplotlib.pyplot as plt

class Particle:
    """
        Particle
    """

    def __init__(self, pos0, vel0, mass=1):
        # current position and velocity
        self.pos = pos0
        self.vel = vel0
        self.mass = mass

        # list for saving all positions of particle
        self.trajectory = []
        self.trajectory.append(self.pos)

        # list for saving all positions also outside the box of particle
        self.trajectory_out = []
        self.trajectory_out.append(self.pos)
        
        # lists for saving all kinetic and total energies
        self.kin_energy = []
        self.tot_energy = []
        


class MolecularDynamics:
    """
    Whole simulation of n particles in a box interacting via Lennard-Jones-potential.
    """

    def __init__(self, rLength, num_particles, time_step, dimension=2, potential='LJ', method='Verlet', new_mic=True):
        self.length = rLength
        self.num_particles = num_particles
        self.h = time_step
        self.dimension = dimension
        self.method = method
        self.new_mic = new_mic

        if potential == 'LJ':
            self.force = forces.LJ_force
            self.potential = forces.LJ_potential_energy
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
        positions = self.length/4 + self.length*0.5 * np.random.rand(self.num_particles, self.dimension)
        velocities = 2 * self.length/10 * np.random.rand(self.num_particles, self.dimension) - self.length/10

        # positions = np.array([[4,6,5],[6,6,5]])
        # velocities = np.array([[0.5,0,0.1],[-0.5,0,-0.1]])
        for i in range(self.num_particles):
            (self.Particles).append(Particle(positions[i,:], velocities[i,:]))
            print("pos: {} \t vel: {}".format(self.Particles[i].pos, self.Particles[i].vel))

        return 0

    def minimum_image_convention(self):
        """
        Minimum image convention for all particles.

        :return: List of list of difference vectors of the (mirrored) particles.
        """
        if self.new_mic:
            diff_vectors = []
            for particle in self.Particles:
                other_particles = self.Particles.copy()
                other_particles.remove(particle)
                loc_diff_vectors = []  # diff vectors for this one particular particle
                for other_particle in other_particles:
                    diff_vector = (particle.pos - other_particle.pos + self.length/2) % self.length - self.length/2
                    loc_diff_vectors.append(diff_vector)
                diff_vectors.append(loc_diff_vectors)
        else:
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

        if self.method == 'Verlet_x':
            if len(self.Particles[0].trajectory) >= 2:
                # normal simulation step
                for i, particle in enumerate(self.Particles):
                    diff_vectors = list_diff_vectors[i]

                    # ------- saving energies first
                    particle.kin_energy.append(forces.kinetic_energy(particle.vel))
                    particle.tot_energy.append(forces.kinetic_energy(particle.vel) + self.potential(diff_vectors))
                    # -------

                    if (np.linalg.norm(particle.pos - particle.trajectory[-2])>self.length/2):
                        print("Alert!")
                        print("using: ", particle.trajectory_out[-1])
                        print("other than: ", particle.trajectory[-1])
                        pos = particle.trajectory_out[-1]
                    else:
                        pos = particle.trajectory[-1]

                    new_position = 2*pos - particle.trajectory[-2] + self.force(diff_vectors) * self.h**2
                    particle.vel = (new_position - particle.trajectory[-2]) /(2*self.h)

                    (particle.trajectory_out).append(new_position)
                    new_position = new_position % self.length

                    (particle.trajectory).append(new_position)
                    particle.pos = new_position
            else:
                # previous timestep unknown, aproximate it
                for i, particle in enumerate(self.Particles):
                    diff_vectors = list_diff_vectors[i]

                    # ------- saving energies first
                    particle.kin_energy.append(forces.kinetic_energy(particle.vel))
                    particle.tot_energy.append(forces.kinetic_energy(particle.vel) + self.potential(diff_vectors))
                    # -------

                    previous_position = particle.pos - self.h * particle.vel
                    (particle.trajectory_out).append(new_position)

                    new_position = 2 * particle.pos - previous_position + self.force(diff_vectors) * self.h ** 2
                    particle.vel = (new_position - previous_position) / (2 * self.h)

                    new_position = new_position % self.length

                    (particle.trajectory).append(new_position)
                    particle.pos = new_position

        elif self.method == 'Verlet_v':
            list_diff_vectors_t = list_diff_vectors
            for i, particle in enumerate(self.Particles):
                diff_vectors = list_diff_vectors_t[i]

                # ------- saving energies first
                particle.kin_energy.append(forces.kinetic_energy(particle.vel))
                particle.tot_energy.append(forces.kinetic_energy(particle.vel) + self.potential(diff_vectors))
                # -------

                new_position = particle.pos + self.h * particle.vel + (self.h ** 2 / 2) * self.force(diff_vectors)
                new_position = new_position % self.length
                particle.trajectory.append(new_position)
                particle.pos = new_position

            list_diff_vectors_th = self.minimum_image_convention()
            for i, particle in enumerate(self.Particles):
                diff_vectors_t = list_diff_vectors_t[i]
                diff_vectors_th = list_diff_vectors_th[i]

                particle.vel += (self.h / 2) * (self.force(diff_vectors_th) + self.force(diff_vectors_t))

        elif self.method == 'Euler':
            for i, particle in enumerate(self.Particles):
                diff_vectors = list_diff_vectors[i]

                # ------- saving energies first
                particle.kin_energy.append(forces.kinetic_energy(particle.vel))
                particle.tot_energy.append(forces.kinetic_energy(particle.vel) + self.potential(diff_vectors))
                # -------

                new_position = particle.pos + particle.vel * self.h
                new_position = new_position % self.length
                particle.vel += self.h * self.force(diff_vectors)

                (particle.trajectory).append(new_position)
                particle.pos = new_position

        else:
            print("Method for simulation step is not implemented.")
            exit(-1)


        return 0

    def simulate(self, num_steps, save_filename=None, save_filename_energies=None):

        for t in range(num_steps):
            self.simulation_step()

        if save_filename is None:
            self.save_filename = "Trajectories.txt"
        else:
            self.save_filename = save_filename
            
        # ----- the same for filename of energies
        if save_filename_energies is None:
            self.save_filename_energies = "Energies.txt"
        else:
            self.save_filename_energies = save_filename_energies
        
        self.energies = self.save_energies(self.save_filename_energies)
        # -----

        self.save_trajectories(self.save_filename)
        return 0

    def save_trajectories(self, save_filename):
        
        trajectories = []
        for particle in self.Particles:
            trajectories.append(particle.trajectory)
        plotting.save_to_file(save_filename, trajectories, self.length)
        
        return trajectories
    
    def save_energies(self, save_filename):
        
        energies = []
        for particle in self.Particles:
            energies.append([particle.kin_energy,particle.tot_energy])
        plotting.save_energies_to_file(save_filename, energies)
        
        return energies

    def plot_trajectories(self):
        trajectories = self.save_trajectories(self.save_filename)
        if self.dimension == 1:
            plotting.plot_1D(trajectories, self.length)
        elif self.dimension == 2:
            plotting.plot_2D(trajectories, self.length)
        elif self.dimension == 3:
            plotting.plot_3D(trajectories, self.length)
        else:
            print("Dimension bigger or equal than 4. No reasonable method for plotting implemented!")

        return 0
    
    def plot_energies(self):  #FINISH ME!!!
        x = np.arange(len(self.energies[0][0]))
        fig, axs = plt.subplots(3)

        energy = 0
        for particle in range(self.num_particles):
            axs[0].plot(x, self.energies[particle][0], label='Particle ' + str(particle), linewidth=2.5)
            axs[1].plot(x, self.energies[particle][1], label='Particle ' + str(particle), linewidth=2.5)
            energy += np.array(self.energies[particle][1])

        axs[2].plot(x, energy)

        axs[0].set_title('Kinetic Energies')
        axs[1].set_title('Total Energies')
        axs[2].set_title('Total Energy of the System')

        axs[0].legend(loc='best')
        axs[1].legend(loc='best')

        plt.tight_layout()
        plt.show()
        return 0
