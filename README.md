# ArgonSimulation
This repository includes a module to simulate molecular dynamics and provides a method simulating a system of N atoms and measuring physical properties such as the pair correlation and the pressure of the system.

By simply pulling the repository you should be able to execute the main.py file. This file contains different important parts of the module:
First the main class "MolecularDynamics" is importet via "from MolDyn import MolecularDynamics". 

## Simulation
Since we want to simulate the system and make measurements of physical properties at the same time, we define the time the simulation will be reinitialized (NUM_INIT) and the time the simulation runs (SIM_TIME). Choosing a large NUM_INIT yields a quite unbiased estimate of the ensemble average of the system. A large SIM_TIME results in measurements at very distinct timestamps, which is desirable since measuring an entity at the very same time does not give any better statistics. However, there is always a tradeoff between accuracy and computational expense.

The basis of *every* simulation would be the creation of an object of class MolecularDynamics with specified temperature and density.
After that, one can run a simulation of the system for a specified time t via object.simulate(t). 

Since we want to measure the physical properties, we will call the member function stat_properties(NUM_INIT, SIM_TIME), which on its own carries out multiple measurements of the pair correlation and the pressure. Via enabling direct plotting via "plot=True", it is directly possible to see the outcomes of the simulation, such as the trajectories of the simulations as well as the measurements.

## Plotting given data
The simulation can take quite a while, depending on how many particles are included, how many initializations are done and for how long the system is simulated. Therefore it might be desirable to only plot the data which is already given.

For carrying out these plots, we use the module plotting. The correct functions for plotting are given in the main.py (commented). These functions take should be called with the correct path to the file, as well as (in some cases) some system constants. These can be seen in the file for the trajectories, but here will be just given for you.

## External constant electric field
The model also includes a possibility to add an external constant electric field to the system. This is allowed by setting the optional parameter "const_electric_field" while initializing the model as a vector of desired constant electric field. Charge of the particles is then defaultly set to 1 (in our simulation units), but can be changed via optional parameter "charge".

The model initializes into equilibrium without the electric force and than "switches on" the electric field only for the simulation itself.

In this case the interesting properties to observe are the trajectories of particles and the change in total energy of the system. In main.py we present the code which simulates the system once with given electric field (EL_FILD) for the time (SIM_TIME2) and plots the trajectories and energies.

Carrying out these plots is again possible externaly as well in the same way as described above.