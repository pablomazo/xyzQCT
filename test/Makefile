FC=gfortran
FFLAGS=-O3 -fcheck=all

tests=testing_utils.f90 test_full_trajs.f90

test: $(tests)
	$(FC) $(FFLAGS) -o test $^

run:
	./test

clean:
	rm test
