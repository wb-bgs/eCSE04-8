World Magnetic Anomaly Model (WMAM)
===================================

The purpose of this code is for fitting global geomagnetic models;
it makes use of the [SLATEC](../../libs/slatec/README.md) library.

WMAM contains many functions and subroutines that were formerly part
of the [GlobLibI](./src/globlibi) library. A few of these routines
can be offloaded to GPU via the setting of preprocessor definitions.
It was found that code profilers such as NVIDIA Nsight Systems do not
see this offloaded code if it exists within a separate GlobLibI library
that is statically linked to WMAM. Hence, the GlobLibI source was
moved directly to the WMAM source folder.

The WMAM code takes as input, WDMAM anomaly data and infers a magnetic field in the form
of a spherical harmonic model. The code implements the conjugate gradient iterative
inversion scheme, one that incorporates the Polak-Ribière formula.

The [build-archer2.sh](build-archer2.sh) script builds WMAM for [ARCHER2](https://www.archer2.ac.uk/)

Simply run `./build-archer2.sh <CPE release> <compiler environment> <build type>` to build the WMAM code.

As of July 2023, the Cray Programming Environment (CPE) release installed on ARCHER2 is `22.12`.
The CPE release version may change in the future - this can be checked by running `module avail cpe`.
If the chosen CPE release is not installed on ARCHER2, the default CPE release will be used instead.

The compiler environment can be `cray`, `gnu` or `aocc`; these are the three compiler
environments provided by the CPE release.

The build type can be `release`, `debug`, `craypat`, `armmap` or `scalasca`; the last three
build types refer to profiling tools, [Cray PAT](https://docs.archer2.ac.uk/user-guide/profile/#craypat), \
[Arm MAP](https://docs.archer2.ac.uk/data-tools/arm-forge/) and [Scalasca](https://www.archer2.ac.uk/training/courses/230822-scalasca/).

Once the build has completed the executable file will be written to
`${HOME/home/work}/apps/wmam/`.
