World Magnetic Anomaly Model (WMAM)
===================================

This library contains general code for fitting global geomagnetic models;
it makes use of [globlibi](../../libs/globlibi/README.md) and the [SLATEC library](../../libs/slatec/README.md).

This model code takes as input WDMAM anomaly data and infers a magnetic field
represented in terms of spherical harmonics. The code offers a choice of two
iterative inversion schemes, Polak-Ribière and conjugate gradient.

The [build.sh](build.sh) script builds WMAM for [ARCHER2](https://www.archer2.ac.uk/)

Simply run `./build.sh <CPE release> <compiler environment> <build type>` to build the WMAM application.

As of Dec 2022, the Cray Programming Environment (CPE) release can be `21.04`, `21.09` or `22.04`.
The available releases may change in the future - this can be checked by running `module avail cpe`.
If the chosen CPE release is not installed on ARCHER2, the default CPE release will be used instead.

The compiler environment can be `cray`, `gnu` or `aocc`; these are the three compiler
environments provided by the CPE release.

The build type can be `release`, `debug`, `craypat`, `armmap` or `scalasca`; the last three
build types refer to profiling tools, [Cray PAT](https://docs.archer2.ac.uk/user-guide/profile/#craypat), [Arm MAP](https://docs.archer2.ac.uk/data-tools/arm-forge/) and [Scalasca](https://www.archer2.ac.uk/training/courses/220509-scalasca/).

Once the build has completed the executable file will be written to
`${HOME/home/work}/apps/wmam/`.
