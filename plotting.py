import numpy as np
import matplotlib.pyplot as plt

def load_from_file_(name_of_file="Trajectories.txt"):
    """
    Loads particle trajectories from a file

    :return: List of list of trajectories of particles and length of the box.
    """
    try:
        file = open(name_of_file,"r")
    except IOError:
        print ("Could not open file.")
    parameters = file.readline() #reads the first line, there is number of particles, number of steps of simulation, dimension of simulation
    parameters = parameters.rstrip('\n')
    parameters = parameters.split(" ")
    number_of_particles = int(parameters[0])
    number_of_steps = int(parameters[1])
    dimension = int(parameters[2])
    length = float(parameters[3])
    all_trajectories=[]

    for particle in range(number_of_particles):
        trajectory = []
        for step in range(number_of_steps):
            position = file.readline()
            position = position.rstrip(' \n')
            position = position.split(" ")
            position = np.array(position, dtype = float)
            trajectory.append(position)
        all_trajectories.append(trajectory)
    file.close()
    return all_trajectories, dimension, length

def save_to_file(name_of_file, trajectories, length):
    """
    Saves particle trajectories from the simulation to a file.
    return: 0 if succes
    """
    number_of_particles = len(trajectories)
    number_of_steps = len(trajectories[0])
    dimension = len(trajectories[0][0])
    file = open(name_of_file,"w")
    file.write("%d %d %d %d\n" %(number_of_particles,number_of_steps,dimension, length))
    for particle in range(number_of_particles):
        for step in range(number_of_steps):
            for i in range(dimension):
                file.write("%f " %trajectories[particle][step][i])
            file.write("\n")
    file.close()
    return 0

def save_energies_to_file(name_of_file, energies):
    """
    Saves particle energies from the simulation to a file.
    return: 0 if success
    """
    number_of_particles = len(energies)
    number_of_steps = len(energies[0][0])
    file = open(name_of_file,"w")
    file.write("%d %d\n" %(number_of_particles,number_of_steps))
    for particle in range(number_of_particles):
        for step in range(number_of_steps):              #saving kinetic energies
            file.write("%f " %energies[particle][0][step])
        file.write("\n")
        for step in range(number_of_steps):              #saving total energies
            file.write("%f " %energies[particle][1][step])
        file.write("\n")
            
    file.close()
    return 0

def load_energies_from_file_(name_of_file="Energies.txt"):
    """
    Loads particle energies from a file.

    :return: List of list of kinetic(0) and total(1) energies of particles.
    """
    try:
        file = open(name_of_file,"r")
    except IOError:
        print ("Could not open file.")
    parameters = file.readline() #reads the first line, there is number of particles, number of steps of simulation, dimension of simulation
    parameters = parameters.rstrip('\n')
    parameters = parameters.split(" ")
    number_of_particles = int(parameters[0])
    #number_of_steps = int(parameters[1])
    all_energies=[]

    for particle in range(number_of_particles):
        energies = []
        kin_en = file.readline()
        kin_en = kin_en.rstrip(' \n')
        kin_en = kin_en.split(" ")
        kin_en = np.array(kin_en, dtype = float)
        energies.append(kin_en)
            
        tot_en = file.readline()
        tot_en = tot_en.rstrip(' \n')
        tot_en = tot_en.split(" ")
        tot_en = np.array(tot_en, dtype = float)
        energies.append(tot_en)
        all_energies.append(energies)
    file.close()
    return all_energies

def plot_2D(trajectories, length):
    color='r'
    for trajectorie in trajectories:
        x = []
        y = []
        for i in range(len(trajectorie)):
            x.append(trajectorie[i][0])
            y.append(trajectorie[i][1])

        split_indices = [0]
        for i in range(len(x)-1):
            if np.abs(x[i] - x[i+1]) >= 0.8*length or \
                    np.abs(y[i] - y[i + 1]) >= 0.8*length:
                split_indices.append(i+1)
        split_indices.append(len(x))

        plot_x = []
        plot_y = []
        for i in range(len(split_indices)-1):
            plot_x.append(x[split_indices[i]:split_indices[i+1]])
            plot_y.append(y[split_indices[i]:split_indices[i+1]])

        for x_part, y_part in zip(plot_x, plot_y):
            plt.plot(x_part, y_part, linewidth=1.2, color=color)

        plt.plot(x[0], y[0], marker='o', markersize=2.0, color='k')

    plt.xlim(0, length)
    plt.ylim(0, length)

    plt.title('Trajectory')
    plt.show()

    return 0

def plot_3D(trajectories, length):
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    color='r'

    for trajectorie in trajectories:
        x = []
        y = []
        z = []
        for i in range(len(trajectorie)):
            x.append(trajectorie[i][0])
            y.append(trajectorie[i][1])
            z.append(trajectorie[i][2])

        split_indices = [0]
        for i in range(len(x) - 1):
            if np.abs(x[i] - x[i + 1]) >= 0.8 * length or \
                np.abs(y[i] - y[i + 1]) >= 0.8 * length or \
                np.abs(z[i] - z[i + 1]) >= 0.8 * length:
                split_indices.append(i + 1)
        split_indices.append(len(x))

        plot_x = []
        plot_y = []
        plot_z = []
        for i in range(len(split_indices) - 1):
            plot_x.append(x[split_indices[i]:split_indices[i + 1]])
            plot_y.append(y[split_indices[i]:split_indices[i + 1]])
            plot_z.append(z[split_indices[i]:split_indices[i + 1]])

        for x_part, y_part, z_part in zip(plot_x, plot_y, plot_z):
            ax.plot(x_part, y_part, z_part, linewidth=1.2, color=color)

        ax.plot(x[0], y[0], z[0], marker='o', markersize=2.0, color='k')
    ax.set_xlim3d(0,length)
    ax.set_ylim3d(0,length)
    ax.set_zlim3d(0,length)

    ax.set_title('Trajectory')
    plt.show()

    return 0

def plot_trajectories(name_of_file):
    trajectories, dimension, length = load_from_file_(name_of_file)
    if dimension == 1:
        plot_1D(trajectories, length)
    elif dimension == 2:
        plot_2D(trajectories, length)
    elif dimension == 3:
        plot_3D(trajectories, length)
    else:
        print("Dimension bigger or equal than 4. No reasonable method for plotting implemented!")

    return 

def plot_energies(name_of_file):  #FINISH ME!!!
    energies = load_energies_from_file_(name_of_file)
    number_of_particles = len(energies)
    number_of_steps = len(energies[0][0])
    print(number_of_particles,number_of_steps)
    x = np.arange(number_of_steps)
    fig, axs = plt.subplots(2)

    for particle in range(number_of_particles):
        axs[0].plot(x, energies[particle][0], label='Particle ' + str(particle),
                    linewidth=2.5)
        axs[1].plot(x, energies[particle][1], label='Particle ' + str(particle),
                    linewidth=2.5)
    axs[0].set_title('Kinetic Energies')
    axs[1].set_title('Total Energies')

    axs[0].legend(loc='best')
    axs[1].legend(loc='best')

    plt.tight_layout()
    plt.show()
    return 0


    