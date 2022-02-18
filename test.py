from MolDyn import MolecularDynamics
import plotting as plotting

model = MolecularDynamics(10, 2, 0.01, dimension=3)

model.simulate(700)
model.plot_trajectories()

#plotting.plot_trajectories("Trajectories.txt")


