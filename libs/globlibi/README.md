Iterative Globel Model Fitting Library (GlobLibI)
=================================================

This library contains general code for fitting global geomagnetic models;
it makes use of the [SLATEC library](../slatec/README.md).

The GlobLibI library is required by the World Magnetic Anomaly Model (WMAM) code.
The [build-archer2.sh](build-archer2.sh) script builds globlibi for [ARCHER2](https://www.archer2.ac.uk/).

Simply run `./build-archer2.sh <CPE release> <compiler environment> <build type>` to build the GlobLibI library.

As of July 2023, the Cray Programming Environment (CPE) release installed on ARCHER2 is `22.12`.
The available CPE release version may change in the future - this can be checked by running `module avail cpe`.
If the chosen CPE release is not installed on ARCHER2, the default CPE release will be used instead.

The compiler environment can be `cray`, `gnu` or `aocc`; these are the three compiler
environments provided by the CPE release.

The build type can be `release`, `debug`, `craypat`, `armmap` or `scalasca`; the last three
build types refer to profiling tools, [Cray PAT](https://docs.archer2.ac.uk/user-guide/profile/#craypat), \
[Arm MAP](https://docs.archer2.ac.uk/data-tools/arm-forge/) and [Scalasca](https://www.archer2.ac.uk/training/courses/230822-scalasca/).

Once the build has completed the library file(s) will be written to
`${HOME/home/work}/libs/globlibi/`.
