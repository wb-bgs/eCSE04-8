# globlibi - an iterative version of the global model fitting library

## Requirements / preparation

This library is only expected to compile on standard Geomag-Team
Linux Virtual Machines (VMs) or the BGS HPC cluster (`bgskwcluster1`).
It is not expected to compile with Solaris.

This library requries a compiler (it's been tested with GNU and
Intel compilers) and appropriate MPI libraries.

When linking with this library to form an executable, the SLATEC
library will also have to be linked as the globlibi library requires
routines from it.

On a standard Geomag-Team Linux VM, you should load the OpenMPI
module with the command

```bash
module load mpi/openmpi-x86_64
```

or equivalent (use `module av` to see what modules are available).

On the BGS HPC cluster load the following libraries with commands

```bash
module load shared \
            slurm \
            dot \
            intel-cluster-runtime/intel64/3.8 \
            intel/compiler/64/2017/17.0.1 \
            intel/mpi/64/2017/1.132
```

or equivalent (again use `module av` to see what's available)
although they may not all be required for compilation.

## Compiling and installing the library

If you wish to check what the make commands are going to do
before doing it (for example to find out where the `make install`
command will put the library), add the `-n` flag, which just prints
commands and does not execute them (e.g. use `make -n install`).

This library can be compiled by simply executing

```bash
make
```

This should place `libgloblibi.a` (the trailing `i` is for
'iterative') at

```bash
./${NERCARCH}/libgloblib.a
```

`NERCARCH` is an environment variable that most likely evaluates
to `linux`.

And if you wish to install, execute

```bash
make install
```

This will place the `libgloblib.a` in a standard location.  Use
the `make -n install` command to see where `make` is planning to
put it.

## Using the library

This library contains general code for fitting global geomagnetic
models but routines will have to be compiled and linked with the
library in order for it to be used.  For example, the files
(not included here)
```bash
build_damp.f
build_damp_space.f
sph_wmam.f
```

add code for the specific problem you want to solve and a driver
program, for example (not included here)

```
mod_wmam_020.f
```

is required to link everything together into an executable.
These files (or others like them) are held separately with the
individual model runs.

You will also have to link with the SLATEC library to produce
executables such as `mod_wmam_020`.
