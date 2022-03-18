from MolDyn import MolecularDynamics
import plotting as plotting

model = MolecularDynamics(1.2, temperature=0.5*119.8, num_unit_cells=2, dimension=3, method='Verlet_v', new_mic=True)

model.stat_properties(1, 3, plot=True)

#plotting.plot_PC("PC_rho=12_T=05_N=32.txt", 2.987603164371443, 3, 32)
#plotting.plot_PC("PC_rho=08_T=10_N=32.txt", 3.4199518933533937, 3, 32)
#plotting.plot_PC("PC_rho=03_T=30_N=32.txt", 4.74252440598675,3,32)
#plotting.plot_PC("PC_rho=08_T=10_N=108.txt", 5.12992784003009, 3, 108)


plotting.plot_trajectories("Trajectories.txt")
# model.simulate(0.1)
#model.plot_trajectories()
#model.plot_energies()
# model.plot_velocity()

#plotting.plot_trajectories("Trajectories.txt")
# plotting.plot_energies("Energies.txt")