FC=gfortran
FFLAGS=-O3

DDEABM=./roots-fortran/src/root_module.F90 ./ddeabm/src/ddeabm_module.F90
DDEABMo=$(DDEABM:%.F90=%.o)
	
# Make ddeabm
ddeabm: $(DDEABMo)
	rm -f libddeabm.a
	ar cvr libddeabm.a $^

clean:
	rm -f *a *mod *.o ./roots-fortran/src/*.o ./ddeabm/src/*.o ./roots-fortran/src/*.mod ./ddeabm/src/*.mod

%.o: %.F90
	$(FC) $(FFLAGS) -o $@ -c $<
