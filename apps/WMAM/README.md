World Magnetic Anomaly Model (WMAM)
===================================

This library contains general code for fitting global geomagnetic models;
it makes use of [globlibi](../../libs/globlibi/README.md) and the [SLATEC library](../../libs/slatec/README.md).

This model code takes as input WDMAM anomaly data and infers a magnetic field
represented in terms of spherical harmonics. The code offers a choice of two
iterative inversion schemes, Polak-Ribière and conjugate gradient.

The [build.sh](build.sh) script builds WMAM for [ARCHER2](https://www.archer2.ac.uk/).
