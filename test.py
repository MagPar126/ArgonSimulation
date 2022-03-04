from MolDyn import MolecularDynamics
import plotting as plotting

model = MolecularDynamics(10, 2, 0.01, dimension=3, method='Verlet_v', new_mic=True)

model.simulate(5000)
model.plot_trajectories()
model.plot_energies()

# plotting.plot_trajectories("Trajectories.txt")
# plotting.plot_energies("Energies.txt")