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

The [build-cirrus-nvfortran.sh](build-cirrus-nvfortran.sh) script builds WMAM for use on
the [Cirrus](https://www.cirrus.ac.uk/) NVIDIA V100 GPUs.

Simply run `./build-cirrus-nvfortran.sh <release|debug>` to build the WMAM code.

Once the build has completed the executable file will be written to
`${HOME/home/work}/apps/wmam/`.
