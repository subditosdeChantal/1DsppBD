# One-dimensional off-lattice system of self-propelled particles.

1. Particles follow Brownian Dynamics with self propulsion.

1. Particle interactions are given by a Lennard-Jones potential.

1. The system has Periodic Boundary Conditions.

1. Particles can also have explicit collective dynamics: moving clusters.

**IMPORTANTE**: Al menos de momento, las dos carpetas del repositorio deben estar a distancia __../..__ de home para que todo funcione correctamente. El repositorio debe colocarse en home por lo tanto.

**INSTALLATION**: For everything to work fine you need to clone this repository into your home directory (__~/__ in Linux). Once this is done you're ready to go. 

**SIMULATION LAUNCHING**: To launch a simulation you just need to execute __./BD__ in __~/1DsppBD/codes/__. This will start a new simulation with the parameters provided in __~/1DsppBD/codes/input.dat__ and create a new directory in __~/output_1DsppBD/__ where it will store all the results of the simulation. To modify the parameters of the simulation you need just edit the __~/1DsppBD/codes/input.dat__ file following the guidelines given in __~/1DsppBD/codes/input_HOWTO.dat__.

**RESULT ANALYSIS AND FIGURE CREATION**: Once you have completed all the simulations you wished, you can just type __bash run_analysis.sh__ in __~/1DsppBD/analysis/__ to process the data of all your simulations. This will create figures for a number of relevant variables of the system for each simulation in their respective output directories.

**LAUNCH WHOLE THING TOGETHER**: You can also set up a simulation loop with __~/1DsppBD/codes/input.dat__ and then run both the simulation and the later analysis codes with a single bash script by typing __bash run.sh__ in __~/1DsppBD/__.
