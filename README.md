# PRISM

This repository contains the source code for the Primitive Equations In Spectral coordinates Model (PRISM) tidal and planetary wave model. The core Fortran77 scripts are organized in the source directory, while useful Python scripts for input and output processing are located in the pyscripts directory. 

PRISM has been used primarily for the simulation of upper atmospheric global-scale waves. Key features are its 3-dimensional background atmospheric nudging capabilities, and its flexible specification of thermal and gravitational tidal forcing terms. 

The current version of the model has been tested with the gfortran compiler on linux ubuntu/hpc systems, but past versions have also been compiled successfully using the PGI compiler. Compiler flags and commands specific to gfortran are included in the source/runmod.sh file. Running the auxiliary .py scripts requires any version of Python 3 to be installed, along with the dependencies described in pyscripts/dependencies.txt. 

The PRISM code can be downloaded by running

```bash
git clone git@github.com:WillemVanCaspel/PRISM.git
```

Once the configuration file (source/input.ini) has been set up, the model can be compiled and run using the following commands:

- make clean
- compile
- sbash modrun.sh

A more comprehensive tutorial on how to run and configure the model using simplified (dummy) input files is given in the README.md file in the tutorial directory.  

Relevant literature
-----------------------------------

 van Caspel, W. E., Espy, P. J., Ortland, D. A., & Hibbins, R. E. (2022). The mid- to high-latitude migrating semidiurnal tide: Results from a mechanistic tide model and SuperDARN observations. Journal of Geophysical Research: Atmospheres, 127, e2021JD036007. https://doi.org/10.1029/2021JD036007 
