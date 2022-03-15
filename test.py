from MolDyn import MolecularDynamics
import plotting as plotting

model = MolecularDynamics(1.2, temperature=0.5*119.8, num_unit_cells=1, dimension=3, method='Verlet_v', new_mic=True)

model.stat_properties(2, 5)

# model.simulate(0.1)
#model.plot_trajectories()
#model.plot_energies()
# model.plot_velocity()

# plotting.plot_trajectories("Trajectories.txt")
# plotting.plot_energies("Energies.txt")