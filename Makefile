FC=gfortran
FFLAGS=-O3

ROOT=./# Path to source files of xyzQCT
QCT=$(addprefix $(ROOT)/src/, constants.f90 settings.f90 hamiltonian.f90 physics.f90 initial_conditions.f90 propagator.f90 utils.f90 QCT.f90 )
QCTo=$(QCT:%.f90=%.o)

DDEABM=$(ROOT)/lib/roots-fortran/src/root_module.F90 $(ROOT)/lib/ddeabm/src/ddeabm_module.F90
DDEABMo=$(DDEABM:%.F90=%.o)
FLIB=-L$(ROOT)/lib/
FMOD=-I$(ROOT)/lib/

#PATH to PES
PESdir=./# Path to potential energy surface
PES=$(addprefix $(PESdir), ) # Potential energy surface source files

# QCT code
QCT: $(PES) $(QCTo) $(ROOT)/src/QCT.f90
	$(FC) $(FFLAGS) -o QCT.x $^ -lddeabm $(FLIB) $(FMOD)

# Make ddeabm
ddeabm: $(DDEABMo)
	rm -f $(ROOT)/lib/ddeabm.a
	ar cvr $(ROOT)/lib/ddeabm.a $^

clean:
	rm -f *.mod *.o *.x $(ROOT)/src/*.mod $(ROOT)/src/*.o $(PESdir)/*.o $(PESdir)/*.mod

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $<