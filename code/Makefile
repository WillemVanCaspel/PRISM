#
#

PROG =	prism.exe

###################################################

include Makefile.SRCS

###################################################

# prefered netCDF 4.2.1.1 or later
LIBS = -lnetcdf -lnetcdff
#explicit pathes needed only if nc-config does not work
INCL = -I/global/hds/software/cpu/eb3/netCDF-Fortran/4.4.4-foss-2017a-HDF5-1.8.18/include
LLIB = -L/global/hds/software/cpu/eb3/netCDF-Fortran/4.4.4-foss-2017a-HDF5-1.8.18/lib

# no netcdf module for prism for now; disabled by setting no path
LIBS = 
INCL = 
LLIB = 

# options by nc-config/nf-config utility
INCL = $(shell nc-config --fflags)
LLIB = $(shell nc-config --flibs)
INCL = 
LLIB =

# gfortran works for both f90 and f77
F90 = gfortran-7

# GNU gfortran compiler (version 4.4.3 or later)
F90FLAGS = -fdec-math -fno-automatic -fdec -ffree-line-length-256 -ffixed-line-length-256 -finit-local-zero
###################################################

LDFLAGS = $(F90FLAGS) $(LLIB) -o $(PROG) $(FOBJ) $(INCL) $(LIBS)

.SUFFIXES: $(SUFFIXES)  .f .f90

.f90.o:
	$(F90) $(F90FLAGS) $(INCL) -c $<

.f.o:
	$(F90) $(F90FLAGS) $(INCL) -c $<

all:  $(PROG)

# Include the dependency-list (created by $ makedepf90 *.f* > dependencies)
include dependencies

$(PROG): $(FOBJ)
	 $(F90) $(LDFLAGS)

clean: diskclean

diskclean:
	rm -f $(PROG) *.o *.mod

##########################################################
