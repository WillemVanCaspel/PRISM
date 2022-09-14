# PRISM

This repository contains the source code for the Primitive Equations In Spectral coordinates Model (PRISM) tidal and planetary wave model. The core Fortran77 code is organized in the source directory, while useful Python scripts for input and output processing are located in the pyscripts directory. 

PRISM has been used primarily for the simulation of upper atmospheric global-scale waves. Key features are its 3-dimensional background atmospheric nudging capabilities, and its flexible specification of thermal and gravitational tidal forcing terms. 

The current version of the model has been tested using the gfortran compiler, but past versions have also been compiled successfully using the PGI compiler. Commands specific to running the model using gfortran (e.g., compiler flags) are included in the source/runmod.sh file. Running the auxiliary .py scripts requires any version of Python 3 to be installed, along with the dependencies described in pyscripts/dependencies.txt. 

```bash
git clone git@github.com:WillemVanCaspel/PRISM.git
```

Once the configuration file source/input.ini has been set up by the user, the model can be run using the following set of commands:

- make clean
- compile
- sbash modrun.sh

By default, simulation results will be sent to the output directory. 

