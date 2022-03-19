#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 12:30:22 2022

@author: Magdalena and Johannes
"""

# how often simulation will be reinitialized to take a statistical average
NUM_INIT = 5
# how long will you simulate
SIM_TIME = 5

from MolDyn import MolecularDynamics

print("\n\n")
###### Simulation of the solid phase ########
model_solid = MolecularDynamics(rho=1.2, temperature=0.5*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True)

model_solid.stat_properties(NUM_INIT, SIM_TIME, plot=True)

print("\n\n")
###### Simulation of the liquid phase ########
model_liquid = MolecularDynamics(rho=0.8, temperature=1.0*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True)

model_liquid.stat_properties(NUM_INIT, SIM_TIME, plot=True)

print("\n\n")
###### Simulation of the gaseous phase ########
model_gas = MolecularDynamics(rho=0.3, temperature=3.0*119.8, num_unit_cells=3, dimension=3, method='Verlet_v', new_mic=True)

model_gas.stat_properties(NUM_INIT, SIM_TIME, plot=True)


