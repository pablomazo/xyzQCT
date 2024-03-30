# xyzQCT

xyzQCT is a code devoted to Quassiclasical Trajectory Calculations (QCT). It implements several initial conditions preparation, such as normal modes and adiabatic switching, as well as propagation of trajectories in order to obtain fixed energy and temperature statistics.

## Download
To clone this repository the following command must be executed:
```
git clone --recurse-submodules git@github.com:pablomazo/xyzQCT.git
```

## Requirements:
- Fortran compiler (the code has been tested with gfortran).
- LAPACK library.
- Compile the DDEABM library provided in the code. For that run:
```
make ddeabm
```

## Usage
In order to use xyzQCT the code must be computed together with the potential energy surface (PES) routines.

xyzQCT will interact with the PES through two subroutines: `setpotxyz` and `potxyz`:
- `setpotxyz`: is called without arguments. This routine will be called at the very beginning of the xyzQCT execution and can be used to load any files or data the PES may require.
- `potxyz(pos, ener, der)`: This routine receives a vector (`pos`) of dimension 3 * number of atoms and returns a float (`ener`) and a vector (`der`) of dimension 3 * number of atoms. `pos` contains the atomic coordinates of all the particles in the system is Bohr. `ener` is the energy for that configuration expressed in Hartree. `der` is the derivative of the energy wrt the positions (negative of the force) and it is expressed in Hartree/Borh.

The files related to the potential energy surface must be included in the Makefile under the variables `PESdir` and `PES`. For example, if the PES files are located in the folder `/home/work/PES/system1` and it is composed of the file `system1_pes.f90` then these variables would be set as:
```
PESdir=/home/work/PES
PES=$(addprefix $(PESdir), system1_pes.f90)
```

After setting these variables run `make` and the `QCT.x` executable should be created.
