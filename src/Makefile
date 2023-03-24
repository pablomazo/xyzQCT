FC=gfortran
FFLAGS=-O3

ROOT=./ # Path to source files of xyzQCT
QCT=$(ROOT)/src/constants.f90 $(ROOT)/src/settings.f90 $(ROOT)/src/QCT.f90 
QCTo=$(QCT:%.f90=%.o)

FLIB=-L$(ROOT)/lib/

#PATH to PES
PESdir= # Path to potential energy surface
PES=$(addprefix $(PESdir), ) # Potential energy surface source files

# QCT code
QCT: $(PES) $(QCTo)
	$(FC) $(FFLAGS) -o QCT.x $^ -lddeabm $(FLIB) 

clean:
	rm -f *.mod *.o *.x $(ROOT)/src/*.mod $(ROOT)/src/*.o $(PESdir)/*.o $(PESdir)/*.mod

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $<
