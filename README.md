# PRISM

This repository contains the source code for the Primitive Equations In Sigma-coordinates Model (PRISM) tidal and planetary wave model. The core Fortran77 scripts are organized in the code directory, while useful Python scripts for input and output processing are located in the pyscripts directory. 

PRISM has been used primarily for the simulation of upper atmospheric global-scale waves. Key features are its 3-dimensional background atmospheric nudging capabilities and its flexible specification of thermal and gravitational tidal forcing terms. 

The current version of the model has been tested with the gfortran compiler on linux systems, but past versions have also been compiled successfully using the PGI compiler. Compiler flags and commands specific to gfortran are included in the code/Makefile file. Running the .py scripts requires any version of Python 3 to be installed, along with the dependencies defined inside the scripts.

The PRISM code can be downloaded using the following command

```bash
git clone git@github.com:WillemVanCaspel/PRISM.git
```

The model can be compiled and run using the following commands from the /code directory:

- make clean
- make
- ./run.sh

This will configure the model as specified by the config.inp file in the main directory.

Relevant literature (selection)
-----------------------------------

 van Caspel, W. E., Espy, P., Hibbins, R., Stober, G., Brown, P., Jacobi, C., & Kero, J. (2023). A case study of the solar and lunar semidiurnal tide response to the 2013 sudden stratospheric warming. Journal of Geophysical Research: Space Physics, 128, e2023JA031680. https://doi.org/10.1029/2023JA031680 

 van Caspel, W. E., Espy, P. J., Ortland, D. A., & Hibbins, R. E. (2022). The mid- to high-latitude migrating semidiurnal tide: Results from a mechanistic tide model and SuperDARN observations. Journal of Geophysical Research: Atmospheres, 127, e2021JD036007. https://doi.org/10.1029/2021JD036007 

 Lieberman, R. S., J. France, D. A. Ortland, and S. D. Eckermann. "The Role of Inertial Instability in Cross-Hemispheric Coupling". Journal of the Atmospheric Sciences 78.4 (2021): 1113-1127. https://doi.org/10.1175/JAS-D-20-0119.1

 D. A. Ortland. Daily estimates of the migrating tide and zonal mean temperature in the mesosphere and lower thermosphere derived from saber data. Journal of Geophysical Research: Atmospheres, 122(7):3754–3785, 2017. doi: https://doi.org/10.1002/2016JD025573

 D. A. Ortland and M. J. Alexander. Gravity wave influence on the global structure of the diurnal tide in the mesosphere and lower thermosphere. Journal of Geophysical Research: Space Physics, 111(A10), 2006. https://doi.org/10.1029/2005JA011467

 D. A. Ortland. A study of the global structure of the migrating diurnal tide using generalized hough modes. Journal of the Atmospheric Sciences, 62(8), 2005. https://doi.org/10.1175/JAS3501.1


