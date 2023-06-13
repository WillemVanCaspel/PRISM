# PRISM

This repository contains the source code for the Primitive Equations In Spectral coordinates Model (PRISM) tidal and planetary wave model. The core Fortran77 scripts are organized in the code directory, while useful Python scripts for input and output processing are located in the pyscripts directory. 

PRISM has been used primarily for the simulation of upper atmospheric global-scale waves. Key features are its 3-dimensional background atmospheric nudging capabilities and its flexible specification of thermal and gravitational tidal forcing terms. 

The current version of the model has been tested with the gfortran compiler on linux systems, but past versions have also been compiled successfully using the PGI compiler. Compiler flags and commands specific to gfortran are included in the code/Makefile file. Running the .py scripts requires any version of Python 3 to be installed, along with the dependencies defined inside the scripts.

The PRISM code can be downloaded by running

```bash
git clone git@github.com:WillemVanCaspel/PRISM.git
```

The model can be compiled and run using the following commands from the /code directory:

- make clean
- make
- ./run.sh

This will configure the model as specified by the config.inp file in the main directory.

Relevant literature
-----------------------------------

 van Caspel, W. E., Espy, P. J., Ortland, D. A., & Hibbins, R. E. (2022). The mid- to high-latitude migrating semidiurnal tide: Results from a mechanistic tide model and SuperDARN observations. Journal of Geophysical Research: Atmospheres, 127, e2021JD036007. https://doi.org/10.1029/2021JD036007 
