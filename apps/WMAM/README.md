World Magnetic Anomaly Model (WMAM)
===================================

This library contains general code for fitting global geomagnetic models;
it makes use of [globlibi](../../libs/globlibi/README.md) and the [SLATEC library](../../libs/slatec/README.md).

This model code takes as input WDMAM anomaly data and infers a magnetic field
represented in terms of spherical harmonics. The code offers a choice of two
iterative inversion schemes, Polak-Ribière and conjugate gradient.

The [build.sh](build.sh) script builds WMAM for [ARCHER2](https://www.archer2.ac.uk/)

Simply run `./build.sh <programming environment> <build type>` to build the WMAM application.

The programming environment can be `cray`, `gnu` or `aocc`; these are the three programming
environments provided by ARCHER2

The build type can be `release`, `debug`, `craypat`, `armmap` or `scalasca`; the last three
build types refer to profiling tools, [Cray PAT](https://docs.archer2.ac.uk/user-guide/profile/#craypat), [Arm MAP](https://docs.archer2.ac.uk/data-tools/arm-forge/) and [Scalasca](https://www.archer2.ac.uk/training/courses/220509-scalasca/).

Once the build has completed the executable file will be written to
`${HOME/home/work}/apps/wmam/`.