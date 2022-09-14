# PRISM

This repository contains the source code for the Primitive Equations In Spectral coordinates Model (PRISM) tidal and planetary wave model. The fortran77 source code is organized in the /source directory, while useful Python scripts for input and output processing are located in the /pyscripts directory. 

PRISM has been used primarily for the simulation of upper atmospheric global-scale waves. Key features are its 3-dimensional background atmosphere nudging capabilities, and its flexible specification of thermal and gravitational tidal forcing terms. 

The current version of the model has been tested to run using the gfortran compiler, but past versions have also been compiled using the xxx compiler. Commands specific to running the model using gfortran (e.g., compiler flags) are described in source/runmod.sh file. Running the auxiliary .py scripts requires Python 3+ to be installed along with the list of dependencies described in the pyscripts/dependencies.txt file. 

Once the configuration file source/input.ini has been set up, the model can be run using the following set (in sequence) of commands:

- make clean
- compile
- sbash modrun.sh

By default, simulation results will be sent to the /output directory. 

