"""
@authors: Magdalena and Johannes

File containing the classes which represents the particle and the tool box for simulating molecular dynamics.

"""

import numpy as np

import forces as forces

import plotting as plotting

import matplotlib.pyplot as plt


class Particle:
    """
        Particle class containing properties like position, velocity and mass.
    """

    def __init__(self, pos0, vel0, mass=1):
        """
        Constructor of particle.

        :param pos0: Initial position of particle.
        :param vel0: Initial velocity of particle.
        :param mass: Mass of particle.
        """
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
    Simulation and measurement of N (same mass) particles in a box interacting via Lennard-Jones-potential.
    """

    ######## INTIALIZATION AND RESETTING ###############
    def __init__(self, rho, temperature, time_step=0.001, dimension=3, num_unit_cells=3,  potential='LJ', method='Verlet_v',
                 new_mic=True):
        """
        Constructor of Molecular Dynamics simulation.

        :param rho: Number density of the system.
        :param temperature: Temperature (in Kelvin) of the system.
        :param time_step: (Optional, default 0.001) Time step for the simulation (<<1).
        :param dimension: (Optional, default 3) Dimension of the system.
        :param num_unit_cells: (Optional, default 3 (=108 particles)) Number of unit cells of the FCC lattice,
         which will be simulated.
        :param potential: (Optional, default 'LJ') String, which tells with which potential particles interact.
         Other potentials needs to be implemented before usage.
        :param method: {'Euler', 'Verlet_x', 'Verlet_v'} (default 'Verlet_v') Method for discretization of EOM and
         simulation step.
        :param new_mic: {'True', 'False'} (default = 'True') Boolean to distinguish an old version of the minimum image
         convention, which was also kind of nice;)
        """
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
        Initialize N particles at the fcc positions with random maxwellian velocities.
        Drive the system to the equilibrium via rescaling the velocities for several timestamps.

        :return: 0 if success
        """

        self.Particles = []

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

            for i in range(num_step):
                self.simulation_step()

            tot_kin = 0
            for particle in self.Particles:
                tot_kin += particle.kin_energy[-1]

            scale_factor = np.sqrt((self.num_particles - 1) * 3 * self.temperature / (2*tot_kin * 119.8))
            print(scale_factor)
            if np.abs(scale_factor - 1) > criterion:
                j +=1
                for particle in self.Particles:
                    particle.vel *= scale_factor
            else:
                for particle in self.Particles:
                    particle.trajectory, particle.trajectory_out, particle.tot_energy, particle.tot_vel, particle.kin_energy =\
                                [particle.trajectory[-1]], [particle.trajectory_out[-1]], [], [particle.tot_vel[-1]], []
                break

        print("Found it! It took {} recursions.".format(j))

        return 0

    def reset(self):
        
        """
        Resets the whole simulation set-up.
        
        :return: 0 if success
        """
        if self.potential == forces.LJ_potential_energy:
            potential = 'LJ'
        else:
            print("Wrong wording or force not implemented yet!")
            exit(-1)
        self.__init__(self.rho, self.temperature, self.h, self.dimension, self.num_unit_cells, potential, self.method,
                      self.new_mic)
        return 0

    def random_initial_conditions(self):
        
        """
        Returns completely random initial conditions (not following physical assumptions).
        
        :return: Initial positions and velocities.
        """
        np.random.seed(0)
        positions = self.length / 4 + self.length * 0.5 * np.random.rand(self.num_particles, self.dimension)
        velocities = 2 * self.length / 10 * np.random.rand(self.num_particles, self.dimension) - self.length / 10
        return positions, velocities

    def initial_positions_face_centered(self):
        
        """
        Sets size of the box based on number density and number of particles.
        Then gives initial positions of atoms in face-centered crystal pattern.
        
        :return: Initial positions of N particles.
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

        return positions

    def initial_maxwellian_velocities(self):
        
        """
        Draws initial random maxwellian velocities for a given temperature of the system.
        
        :return: Non-scaled random maxwellian velocities.
        """
        velocities = np.zeros((self.num_particles,self.dimension))

        velocities = np.random.normal( scale = 2 * self.temperature/119.8, # for argon
            size=(self.num_particles, self.dimension))
    
        from scipy.stats import skew
 
        total_vel_x = sum(velocities[:,0])
        total_vel_y = sum(velocities[:,1])
        total_vel_z = sum(velocities[:,2])
        
        for i,v in enumerate(velocities):
            velocities[i,0] = v[0]  - total_vel_x/self.num_particles
            velocities[i,1]  = v[1] - total_vel_y/self.num_particles
            velocities[i,2]  = v[2] - total_vel_z/self.num_particles
            
        print("means: \t {0:.3f}\t {1:.3f}\t {2:.3f}".format(np.mean(velocities[:,0]), np.mean(velocities[:,1]), np.mean(velocities[:,2])))
        print("skewness: {0:.3f}\t {1:.3f}\t {2:.3f}".format(skew(velocities[:,0]), skew(velocities[:,1]), skew(velocities[:,2])))
        return velocities

    ########## SIMULATION ###############
    def minimum_image_convention(self):
        
        """
        Minimum image convention for all particles.
        Version depends on whether new_mic = True/False. 'True' should always be preferred.

        :return: List of lists of difference vectors of the (mirrored) particles.
        """
        if self.new_mic:
                
            diff_vectors = diff_vectors = [[[] for i in range(self.num_particles-1)] for j in range(self.num_particles)]
            for i,particle in enumerate(self.Particles):
                other_particles = self.Particles.copy()
                other_particles.remove(particle)
                for j in range(i,self.num_particles-1):
                    diff_vector = (particle.pos - other_particles[j].pos + self.length / 2) % self.length - self.length / 2
                    diff_vectors[i][j] = diff_vector
                    diff_vectors[j+1][i] = -diff_vector
            
            
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
        
        """
        Realization of one simulation step according to the given method.
            'Euler' is a standard Euler algorithm,
            'Verlet_x' is a Verlet algorithm in spatial domain and
            'Verlet_v' is a method for Velocity Verlet algorithm.
        Saves energies of all particles, updates their current positions and velocities and saves old positions to trajectories.
        
        :return: 0 if success
        """
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
                # previous timestep unknown, approximate it
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

    def simulate(self, num_time_intervals, save_filename=None, save_filename_energies=None, plot=False, n_means=5):
        
        """
        Executes simulation for a given time.
        Measures pressure and pair correlation n_means times and saves trajectories and energies to a file.

        :param num_time_intervals: Number of time intervals in intrinsic system time.
        :param save_filename: (Optional) Filename to save trajectories of particles. By default, constructed through the
         system parameters.
        :param save_filename_energies: (Optional) Filename to save energies of particles. By default, constructed through
         the system parameters.
        :param plot: (Optional, default = False) Boolean, whether trajectories and energies should be plotted.
        :param n_means: (Optional, default = 5) Number of measurements after half of the simulation time.
        
        :return: 0 if success
        """

        num_steps = int(num_time_intervals/self.h)

        # number measurements per simulation
        t_stamps = np.zeros(n_means, int)
        for i in range(n_means):
            t_stamps[i] = int(num_steps/2 + i*(num_steps/(2*n_means)))
        for t in range(num_steps):
            self.simulation_step()
            if t in t_stamps:
                self.measurement()
                print("Measured.")
            
        print("Done simulating! Now plotting.")

        if save_filename is None:
            self.save_filename = "Trajectories_rho=" + str(self.rho).replace('.', '') + "_T=" + str(self.temperature/119.8).replace('.', '')\
                            + "_N=" + str(self.num_particles) + ".txt"
        else:
            self.save_filename = save_filename

        # ----- the same for filename of energies
        if save_filename_energies is None:
            self.save_filename_energies = "Energies_rho=" + str(self.rho).replace('.', '') + "_T=" + str(self.temperature/119.8).replace('.', '')\
                            + "_N=" + str(self.num_particles) + ".txt"
        else:
            self.save_filename_energies = save_filename_energies

        self.energies = self.save_energies(self.save_filename_energies)
        # -----

        self.save_trajectories(self.save_filename)
        
        if plot == True:
            self.plot_trajectories()
            #self.plot_energies()
            
        return 0

    ########## MEASUREMENT ###############
    def stat_properties(self, num_init, num_time_intervals, save_filename=None, save_filename_energies=None, plot=False):
        
        """
        Measures physical properties (pair correlation and pressure) of a given system.

        :param num_init: Integer which defines how often system and simulation will be reinitialized to take statistical
        average.
        :param num_time_intervals: Number of time intervals for on simulation.
        :param save_filename: (Optional) Filename to save trajectories.
        :param save_filename_energies: (Optional) Filename to sacve energies.
        :param plot: (Optional, default = False) Boolean to decide whether to plot pressure and correlation or not.
        
        :return: 0 if success
        """
        for i in range(num_init):
            print("\n\nSimulation {} out of {}".format(i+1, num_init))
            if i!=0: self.reset()
            self.simulate(num_time_intervals, save_filename, save_filename_energies)
            self.plot_trajectories()


        if plot == True:
            filenamePC = "PC_rho=" + str(self.rho).replace('.', '') + "_T=" + str(self.temperature/119.8).replace('.', '')\
                            + "_N=" + str(self.num_particles) + ".txt"
            plotting.plot_PC(filenamePC, self.length, self.dimension,self.num_particles)
            filenamePP = "PP_rho=" + str(self.rho).replace('.', '') + "_T=" + str(self.temperature / 119.8).replace('.','')\
                            + "_N=" + str(self.num_particles) + ".txt"
            plotting.plot_PP(filenamePP)

    def measure_pressure(self, filenamePP=None):
        
        """
        Measures and saves pressures of the system for one timestamp.

        :param filenamePP: (Optional) Filename to save the pressure.
        
        :return: 0 if success
        """
        
        diff_vectors = self.minimum_image_convention()

        pressure = self.rho*self.temperature/119.8
        for i in range(len(diff_vectors)):
            for j in range(i, len(diff_vectors[i])):
                r = np.linalg.norm(diff_vectors[i][j])
                pressure += (self.rho/(3*self.num_particles)) * (1/2) * 4 * (12*r**(-13) - 6*r**(-7))

        if filenamePP == None:
            filenamePP = "PP_rho=" + str(self.rho).replace('.', '') + "_T=" + str(self.temperature / 119.8).replace('.','')\
                            + "_N=" + str(self.num_particles) + ".txt"

        try:
            file = open(filenamePP, "a")
        except IOError:
            print("Could not open file.")

        file.write(str(pressure))
        file.write("\n")
        file.close()
        return 0

    def measure_corr(self, filenamePC=None):
                
        """
        Measures and saves pair correlation of the system for one timestamp.

        :param filenamePC: (Optional) Filename to save tje correlations.
        
        :return: 0 if success
        """
        diff_vectors = self.minimum_image_convention()
        lengths = []
        for p in diff_vectors:
            for dv in p:
                lengths.append(np.linalg.norm(dv))

        if filenamePC == None:
            filenamePC = "PC_rho=" + str(self.rho).replace('.', '') + "_T=" + str(self.temperature/119.8).replace('.', '')\
                            + "_N=" + str(self.num_particles) + ".txt"

        try:
            file = open(filenamePC, "a")
        except IOError:
            print("Could not open file.")

        file.write(str(lengths))
        file.write("\n")
        file.close()
        return 0

    def measurement(self):
                
        """
        Measures and saves statistical properties of the system.
        
        :return: 0 if success
        """
        self.measure_corr()
        self.measure_pressure()
        return 0

    ########## SAVE RESULTS ##############
    def save_trajectories(self, save_filename):
        
        """
        Saves trajectories of all particles in one big list.

        :param save_filename: Filename to save the trajectories.
        
        :return: Trajectories (List of lists of np.array(dimension)).
        """

        trajectories = []
        for particle in self.Particles:
            trajectories.append(particle.trajectory)
        plotting.save_to_file(save_filename, trajectories, self.length)

        return trajectories

    def save_energies(self, save_filename):
        
        """
        Saves total energies of all particles in one big list.

        :param save_filename: Filename to save the energies.
        
        :return: Total energies (List of floats).
        """

        energies = []
        for particle in self.Particles:
            energies.append([particle.kin_energy, particle.tot_energy])
        plotting.save_energies_to_file(save_filename, energies)

        return energies

    ########### PLOTTING (similar to extra plotting tool, just as methods of the class) #################
    def plot_trajectories(self):
        
        """
        Plots trajectories from a given simulation.

        .. warning:: Only call after simulate().
        
        :return: 0 if success
        """
        trajectories = self.save_trajectories(self.save_filename)
        if self.dimension == 2:
            plotting.plot_2D(trajectories, self.length)
        elif self.dimension == 3:
            plotting.plot_3D(trajectories, self.length)
        else:
            print("Dimension either 1, or bigger than 3. No reasonable method for plotting implemented!")

        return 0

    def plot_energies(self):
            
        """
        Plots total energies from a given simulation.

        .. warning:: Only call after simulate().
        
        :return: 0 if success
        """
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

