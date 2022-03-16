from MolDyn import MolecularDynamics
import plotting as plotting

model = MolecularDynamics(0.3, temperature=3.0*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True)

model.stat_properties(2, 5, plot=True)

# model.simulate(0.1)
#model.plot_trajectories()
#model.plot_energies()
# model.plot_velocity()

# plotting.plot_trajectories("Trajectories.txt")
# plotting.plot_energies("Energies.txt")