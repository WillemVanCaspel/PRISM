# PRISM

This repository contains the source code for the Primitive Equations In Spectral coordinates Model (PRISM) tidal and planetary wave model. The fortran77 source code is organized in the /source directory, while auxiliary Python scripts used for input and output processing are located in the /pyscripts directory. 

PRISM has primarily been used for the simulation of upper atmospheric global-scale waves. Key features are its capability to be nudged to 3-dimensional background meteorological fields, and its flexible specification of thermal and gravitational tidal forcing terms. 

The current version of the model has only been tested to run using the gfortran compiler, but in the past it has also been compiled using .sig. Commands specific to running the model using gfortran (e.g., compiler flags) are located in source/runmod.sh. Running the auxiliary .py scripts requires Python 3+ to be installed along with the list of dependencies described in the pyscripts/dependencies.txt file. 

Once the configuration file source/input.ini has been set up, the model can be run using the following set (in sequence) of commands:

- make clean
- compile
- sbash modrun.sh

By default, simulation results will be sent to the /output directory. 

