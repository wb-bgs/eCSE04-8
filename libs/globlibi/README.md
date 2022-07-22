Iterative Globel Model Fitting Library (globlibi)
=================================================

This library contains general code for fitting global geomagnetic models;
it makes use of the [SLATEC library](../slatec/README.md).

The globlibi library is required by the World Magnetic Anomaly Model code (WMAM).
The [build.sh](build.sh) script builds globlibi for [ARCHER2](https://www.archer2.ac.uk/).

Simply run `./build.sh <programming environment> <build type>` to build the globlibi library.

The programming environment can be `cray`, `gnu` or `aocc`; these are the three programming
environments provided by ARCHER2

The build type can be `release`, `debug`, `craypat`, `armmap` or `scalasca`; the last three
build types refer to profiling tools, [Cray PAT](https://docs.archer2.ac.uk/user-guide/profile/#craypat), [Arm MAP](https://docs.archer2.ac.uk/data-tools/arm-forge/) and [Scalasca](https://www.archer2.ac.uk/training/courses/220509-scalasca/).

Once the build has completed the library file(s) will be written to
`${HOME/home/work}/libs/globlibi/`.