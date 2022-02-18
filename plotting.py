import numpy as np
import matplotlib.pyplot as plt

def load_from_file_(name_of_file):
    """
    Loads particle trajectories from a file

    :return: List of list of trajectories of particles.
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
    #dimension = int(parameters[2])
    all_trajectories=[]

    for particle in range(number_of_particles):
        trajectory = []
        for step in range(number_of_steps):
            position = file.readline()
            position = position.rstrip('\n')
            position = position.split(" ")
            position = np.array(position,dtype = float)
            trajectory.append(position)
        all_trajectories.append(trajectory)
    file.close()
    return all_trajectories

def save_to_file(name_of_file, trajectories):
    """
    Saves particle trajectories from the simulation to a file.
    return: 0 if succes
    """
    number_of_particles = len(trajectories[:,0])
    number_of_steps = len(trajectories[0,:])
    dimension = len(trajectories[0,0])
    file = open(name_of_file,"w")
    file.write("%d %d %d\n" %(number_of_particles,number_of_steps,dimension))
    for particle in range(number_of_particles):
        for step in range(number_of_steps):
            for i in range(dimension):
                file.write("%f " %trajectories[particle,step][i])
            file.write("\n")
    file.close()
    return 0

def plot_2D(trajectories):
    pass

def plot_3D(trajectories):
    pass