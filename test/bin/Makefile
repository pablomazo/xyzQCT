FC=gfortran
FFLAGS=-O3 -fcheck=all

ROOT=../../
QCT=$(addprefix $(ROOT)/src/, xyzqct_constants.f90 xyzqct_utils.f90 xyzqct_lapack.f90 xyzqct_rand.f90 xyzqct_settings.f90 xyzqct_physics.f90 xyzqct_hamiltonian.f90 xyzqct_initial_conditions.f90 xyzqct_propagator.f90)
QCTo=$(QCT:%.f90=%.o)

DDEABM=$(ROOT)/lib/roots-fortran/src/root_module.F90 $(ROOT)/lib/ddeabm/src/ddeabm_module.F90
DDEABMo=$(DDEABM:%.F90=%.o)
FLIB=-L$(ROOT)/lib/ -llapack
FMOD=-I$(ROOT)/lib/

#PATH to PES
PESdir=./
PES=$(addprefix $(PESdir), pot_interface.f90) # Potential energy surface source files

# QCT code
QCT: $(PES) $(QCTo) $(ROOT)/src/qct.f90
	$(FC) $(FFLAGS) -o QCT.x $^ -lddeabm $(FLIB) $(FMOD)

# Make ddeabm
ddeabm: $(DDEABMo)
	rm -f $(ROOT)/lib/ddeabm.a
	ar cvr $(ROOT)/lib/ddeabm.a $^

clean:
	rm -f *.mod *.o *.x $(ROOT)/src/*.mod $(ROOT)/src/*.o $(PESdir)/*.o $(PESdir)/*.mod

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(FMOD)

%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $<
