# hydro1d 

*A Fortran 95 1-d hydrodynamics code*

http://zingale.github.io/hydro1d/


## Problems

Individual problems are built in subdirectories.  For example,
to build and run the Sod shock tube problem, do:

* `cd hydro1d/sod`

* `make`

* `./hydro1d inputs-sod-xp`

As the code is built, object files and modules will be output into the
`_build` subdirectory.  `make realclean` will clean up the objects.
Things are setup for gfortran by default -- you will need to edit the
`Ghydro.mak` with different options for different compilers.  Some bits
of Fortran 2003 and 2008 are used, so an up-to-date compiler is
needed.

## Runtime Parameters

A number of runtime options can be specified -- look in `params.f90`
for the available runtime parameters.  These are set in a namelist in
the inputs file.  Problems can specify their own runtime parameters
(with a `probparams.f90` in the problem directory).  These have a
separate namelist in the same inputs file.  See the sod problem for
examples.

## Development

The PPM implementation should be synced up to what was done in 
Castro's PPM.

