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
        self.tot_vel = []
        self.tot_vel.append(self.vel)


class MolecularDynamics:
    """
    Whole simulation of n particles in a box interacting via Lennard-Jones-potential.
    """

    def __init__(self, rho, time_step=0.001, temperature=100, dimension=3, num_unit_cells=3,  potential='LJ', method='Verlet',
                 new_mic=True):
        self.rho = rho
        self.temperature = temperature
        self.num_unit_cells = num_unit_cells
        self.dimension = dimension
        self.num_particles = int(self.num_unit_cells**self.dimension * (2**self.dimension / 2))
        print("Num particles: ", self.num_particles)
        # system length
        self.length = (self.num_particles/self.rho)**(1/self.dimension)
        print("Lenght: ", self.length)
        # length of the unit cell (as a cube)
        self.unit_size = self.length/(2 * (self.num_unit_cells))
        print("unit size: ", self.unit_size)
        self.h = time_step
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

        #positions, velocities = self.random_initial_conditions()

        positions = self.initial_positions_face_centered()
        velocities = self.initial_maxwellian_velocities()

        sim_time = 0.1 #together with precision 0.05 see
        num_step = int(sim_time/self.h)

        print("Find equilibrium...")
        criterion = 0.1 #seems like 0.05 is quite a reachable precision
        scale_factor = 1.0
        j = 0
        for i in range(self.num_particles):
                (self.Particles).append(Particle(positions[i, :], velocities[i, :]))
        while True:
            #velocities = scale_factor*velocities
            #for i in range(self.num_particles):
                #(self.Particles).append(Particle(positions[i, :], scale_factor * velocities[i, :]))
                #print(self.Particles[i].vel)

            for i in range(num_step):
                self.simulation_step()

            tot_kin = 0
            for particle in self.Particles:
                tot_kin += particle.kin_energy[-1]

            scale_factor *= np.sqrt((self.num_particles - 1) * 3 * self.temperature / (2*tot_kin * 119.8))
            print(scale_factor)
            if np.abs(scale_factor - 1) > criterion:
                j +=1
                for particle in self.Particles:
                    particle.vel *= scale_factor
                
            else:
                break

        print("Found it! It took {} recursions.".format(j))

        return 0

    def random_initial_conditions(self):  # old "dumm" version
        np.random.seed(0)
        positions = self.length / 4 + self.length * 0.5 * np.random.rand(self.num_particles, self.dimension)
        velocities = 2 * self.length / 10 * np.random.rand(self.num_particles, self.dimension) - self.length / 10
        return positions, velocities

    def initial_positions_face_centered(self):
        """
        Sets size of the box based on number density and number of particles . Then gives initial positions of atoms in face-centered crystal patern.
        Returns: initial positions of N particles
        """
        l = self.unit_size
        initial_pos = l/2 * np.ones(self.dimension)
        positions = []

        basis = range(self.num_unit_cells)
        dims = []
        for j in range(self.dimension):
            dims.append(len(basis))

        image_boxes = np.zeros(dims)
        it = np.nditer(image_boxes, flags=['multi_index'])
        for x in it:
            displacement = np.array(
                it.multi_index) * 2 * l  # displacements of copies in respect to the original one
            positions.append(initial_pos + displacement)

        basis = [0,1]
        dims = []
        for j in range(self.dimension):
            dims.append(len(basis))

        image_boxes = np.zeros(dims)
        it = np.nditer(image_boxes, flags=['multi_index'])
        dvs = []
        for x in it:
            displacement = np.array(
                it.multi_index)  # displacements of copies in respect to the original one
            if (np.abs(np.linalg.norm(displacement)**2 - (self.dimension - 1)) < 1e-6):
                dvs.append(l*displacement)

        basis_pos = positions
        filling = [(position + dv ) % self.length for position in positions for dv in dvs]

        positions = np.array(basis_pos + filling)

        # fcc visualization
        # ar = np.array(positions)
        # print(np.amax(ar, axis=0))

        # fig = plt.figure(figsize=(12, 12))
        # ax = fig.add_subplot(projection='3d')
        # ax.set_xlim(0,self.length)
        # ax.set_ylim(0,self.length)
        # ax.set_zlim(0,self.length)

        # ax.scatter(ar[:, 0], ar[:, 1], ar[:, 2])
        # plt.show()

        return positions

    def initial_maxwellian_velocities(self):
        """
        Returns: non-scaled random maxwellian velocities
        """
        #velocities = np.zeros((self.num_particles,self.dimension))
        #for i in range(self.dimension):
        #    vel = np.random.normal( scale = 2 * self.temperature/119.8, # for argon
        #    size=self.num_particles)
        #    for j in range(self.num_particles):
        #        velocities[j,i]=vel[j]
        velocities = np.random.normal( scale = 2 * self.temperature/119.8, # for argon
            size=(self.num_particles, self.dimension))  # Can I do that? Put all direction together?
        return velocities

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
                    diff_vector = (particle.pos - other_particle.pos + self.length / 2) % self.length - self.length / 2
                    loc_diff_vectors.append(diff_vector)
                diff_vectors.append(loc_diff_vectors)
        else:
            diff_vectors = []
            for i, particle in enumerate(self.Particles):
                other_particles = self.Particles.copy()
                other_particles.remove(particle)
                loc_diff_vectors = []  # diff vectors for this one particular particle
                for other_particle in other_particles:
                    copies_distances = []  # list of distances of all copies from the particle
                    copies_diff_vectors = []  # list of diff vectors of all copies from the particle

                    basis = [-1, 0, 1]
                    dims = []
                    for j in range(self.dimension):
                        dims.append(len(basis))

                    image_boxes = np.zeros(dims)
                    it = np.nditer(image_boxes, flags=['multi_index'])
                    for x in it:
                        displacement = (np.array(
                            it.multi_index) - 1) * self.length  # displacements of copies in respect to the original one
                        mirrored_pos = other_particle.pos + displacement  # possition of the mirrored particle
                        diff_vector = particle.pos - mirrored_pos  # diff vector of the mirrored particle
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

                    if (np.linalg.norm(particle.pos - particle.trajectory[-2]) > self.length / 2):
                        print("Alert!")
                        print("using: ", particle.trajectory_out[-1])
                        print("other than: ", particle.trajectory[-1])
                        pos = particle.trajectory_out[-1]
                    else:
                        pos = particle.trajectory[-1]

                    new_position = 2 * pos - particle.trajectory[-2] + self.force(diff_vectors) * self.h ** 2
                    velocity = (new_position - particle.trajectory[-2]) / (2 * self.h)
                    particle.vel = velocity
                    particle.tot_vel.append(velocity)

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

                    new_position = 2 * particle.pos - previous_position + self.force(diff_vectors) * self.h ** 2
                    (particle.trajectory_out).append(new_position)
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

                velocity = particle.vel + (self.h / 2) * (self.force(diff_vectors_th) + self.force(diff_vectors_t))
                particle.vel = velocity
                particle.tot_vel.append(velocity)

        elif self.method == 'Euler':
            for i, particle in enumerate(self.Particles):
                diff_vectors = list_diff_vectors[i]

                # ------- saving energies first
                particle.kin_energy.append(forces.kinetic_energy(particle.vel))
                particle.tot_energy.append(forces.kinetic_energy(particle.vel) + self.potential(diff_vectors))
                # -------

                new_position = particle.pos + particle.vel * self.h
                new_position = new_position % self.length
                velocity = particle.vel + self.h * self.force(diff_vectors)
                particle.vel = velocity

                (particle.trajectory).append(new_position)
                particle.pos = new_position
                (particle.tot_vel).append(velocity)


        else:
            print("Method for simulation step is not implemented.")
            exit(-1)

        return 0

    def simulate(self, num_time_intervals, save_filename=None, save_filename_energies=None):

        num_steps = int(num_time_intervals/self.h)

        for t in range(num_steps):
            self.simulation_step()
            
        print("Done simulating! Now plotting.")

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
            energies.append([particle.kin_energy, particle.tot_energy])
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

    def plot_energies(self):  # FINISH ME!!!
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

    def plot_velocity(self):

        velocities = np.zeros((len(self.Particles[0].tot_vel), 3))
        for particle in self.Particles:
            velocities += particle.tot_vel
        print("Variance of total momentum: ", np.var(velocities, axis=0))
        plotting.plot_3D([velocities], self.length, title='Momentum conservation')
