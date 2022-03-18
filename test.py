from MolDyn import MolecularDynamics
import plotting as plotting

model = MolecularDynamics(0.8, temperature=1.0*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True)

model.stat_properties(5, 5, plot=True)

#plotting.plot_PC("PC_rho=12_T=05_N=32.txt", 2.987603164371443, 3, 32)
#plotting.plot_PC("PC_rho=08_T=10_N=32.txt", 3.4199518933533937, 3, 32)
#plotting.plot_PP("PP_rho=03_T=30_N=108.txt")
#plotting.plot_PC("PC_rho=03_T=30_N=108.txt", 7.113786608980125, 3, 108)


#plotting.plot_trajectories("Trajectories.txt")
# model.simulate(0.1)
#model.plot_trajectories()
#model.plot_energies()
# model.plot_velocity()

plotting.plot_trajectories("Trajectories.txt")
# plotting.plot_energies("Energies.txt")