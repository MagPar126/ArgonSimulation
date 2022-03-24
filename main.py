#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:30:22 2022

@author: Magdalena and Johannes
"""

from MolDyn import MolecularDynamics
import plotting as plotting

# how often simulation will be reinitialized to take a statistical average
NUM_INIT = 5
# how long will you simulate
SIM_TIME = 5


###### Simulation of the solid phase ########
print("\n\nRUNNING: Simulation of the solid phase\n\n")

model_solid = MolecularDynamics(rho=1.2, temperature=0.5*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True)

model_solid.stat_properties(NUM_INIT, SIM_TIME, plot=True)

#plotting.plot_trajectories("Relevant measured data/Trajectories_rho=12_T=05_N=108.txt")
#plotting.plot_PC("Relevant measured data/PC_rho=12_T=05_N=108.txt", 4.481405, 3, 108)
#plotting.plot_PP("Relevant measured data/PP_rho=12_T=05_N=108.txt")


###### Simulation of the liquid phase ########
print("\n\nRUNNING: Simulation of the liquid phase\n\n")

model_liquid = MolecularDynamics(rho=0.8, temperature=1.0*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True)

model_liquid.stat_properties(NUM_INIT, SIM_TIME, plot=True)

#plotting.plot_trajectories("Relevant measured data/Trajectories_rho=08_T=10_N=108.txt")
#plotting.plot_PC("Relevant measured data/PC_rho=08_T=10_N=108.txt", 5.129927, 3, 108)
#plotting.plot_PP("Relevant measured data/PP_rho=08_T=10_N=108.txt")


###### Simulation of the gaseous phase ########
print("\n\nRUNNING: Simulation of the gaseous phase\n\n")

model_gas = MolecularDynamics(rho=0.3, temperature=3.0*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True)

model_gas.stat_properties(NUM_INIT, SIM_TIME, plot=True)

#plotting.plot_trajectories("Relevant measured data/Trajectories_rho=03_T=30_N=108.txt")
#plotting.plot_PC("Relevant measured data/PC_rho=03_T=30_N=108.txt", 7.113787, 3, 108)
#plotting.plot_PP("Relevant measured data/PP_rho=03_T=30_N=108.txt")

#----------------------------------------------------------------------------------------------------------------------#
###### Observing the effects of an external constant electric fiel ########
SIM_TIME2 = 2

print("\n\nRUNNING: Simulation of the solid phase with electric field\n\n")

model_solid_el = MolecularDynamics(rho=1.2, temperature=0.5*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True,const_electric_field = [0.2,0,0])

model_solid_el.simulate(SIM_TIME2, plot=True)

#plotting.plot_trajectories("Relevant measured data/Trajectories_rho=12_T=05_N=108_E=0_2.txt")


print("\n\nRUNNING: Simulation of the solid phase with electric field\n\n")

model_liquid_el = MolecularDynamics(rho=0.8, temperature=1.0*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True,const_electric_field = [2.5,0,0])

model_liquid_el.simulate(SIM_TIME2, plot=True)

#plotting.plot_trajectories("Relevant measured data/Trajectories_rho=12_T=05_N=108_E=2_5.txt")

print("\n\nRUNNING: Simulation of the solid phase with electric field\n\n")

model_gas_el = MolecularDynamics(rho=0.3, temperature=3.0*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True,const_electric_field = [5,0,0])

model_gas_el.simulate(SIM_TIME2, plot=True)

#plotting.plot_trajectories("Relevant measured data/Trajectories_rho=12_T=05_N=108_E=5.txt")





