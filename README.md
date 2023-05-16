# eCSE04-8

This repo holds the software and associated scripts used for the eCSE project,
ARCHER2-eCSE04-8, *Enabling a better global view of the magnetic field of Earth's rocks*.

It contains the (Fortran) source for one code and two libraries.

1. [World Magnetic Anomaly Model (WMAM)](./apps/WMAM/README.md)
2. [Iterative Globel Model Fitting Library (GlobLibI)](./libs/globlibi/README.md)
3. [SLATEC Common Mathematical Library 4.1](./libs/slatec/README.md)

The WMAM code is statically linked to the GlobLibI and SLATEC libraries. This means
you will need to build the libraries before you build WMAM. Build instructions for each
software component can be found via the links above.

Example submission scripts can be found under [./tests/WMAM/submit](./tests/WMAM/submit).
Input data for the *L* = 200, 300 and 720 cases (where *L* is the spherical harmonic degree)
can be found under [./tests/WMAM/shd](./tests/WMAM/shd).

Please note, it is recommended that you first copy the `./tests/WMAM` folder to
`${HOME/home/work}/tests/WMAM`, submitting your jobs from a location outside the
eCSE04-8 repo directory.


This repo also contains utility software.

1. [pypp](./utils/pypp)

   A Python environment for running the WMAM post-processing python scripts,
   see [/tests/WMAM/scripts](/tests/WMAM/scripts).

2. [gmt](./utils/gmt) and [gd2gc](./utils/gd2gc)

   Tools for manipulating geographic and Cartesian data sets. These are required
   for generating the input data for large values of *L*, e.g., *L* = 1440, 2000.
