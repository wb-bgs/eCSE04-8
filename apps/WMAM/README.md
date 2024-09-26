World Magnetic Anomaly Model (WMAM)
===================================

The purpose of this code is for fitting global geomagnetic models;
it makes use of the [GlobLibI](../../libs/globlibi/README.md) and [SLATEC](../../libs/slatec/README.md) libraries.

The WMAM code takes as input, WDMAM anomaly data and infers a magnetic field
in the form of a spherical harmonic model. The code offers a choice of two
flavours of the conjugate gradient iterative inversion scheme; one that
that incorporates the Polak-Ribière formula and one that does not.

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
