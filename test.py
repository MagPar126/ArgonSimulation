#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 23:33:49 2022

@author: majda

"""
from MolDyn import MolecularDynamics
import plotting as plotting
#model_solid = MolecularDynamics(rho=1.2, temperature=0.5*119.8, num_unit_cells=2, dimension=3, method='Verlet_v', new_mic=True,const_electric_field = [2,0,0])

#model_solid.simulate(2,plot=True)
#plotting.plot_trajectories("Trajectories_rho=12_T=05_N=32.txt")
#plotting.plot_trajectories("Relevant measured data/Trajectories_rho=12_T=05_N=108_E=5.txt")
plotting.plot_PC("Relevant measured data/PC_rho=03_T=30_N=108.txt", 7.113787, 3, 108)
plotting.plot_PP("Relevant measured data/PP_rho=03_T=30_N=108.txt")