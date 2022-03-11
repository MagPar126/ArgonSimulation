from MolDyn import MolecularDynamics
import plotting as plotting

model = MolecularDynamics(1.2, temperature=0.5, num_unit_cells=2, dimension=3, method='Verlet_v', new_mic=True)

model.simulate(5)
model.plot_trajectories()
model.plot_energies()
#model.plot_velocity()

# plotting.plot_trajectories("Trajectories.txt")
# plotting.plot_energies("Energies.txt")