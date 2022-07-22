# eCSE04-8

This repo holds the software and associated scripts used for the eCSE project,
ARCHER2-eCSE04-8, *Enabling a better global view of the magnetic field of Earth's rocks*.

It contains the (Fortran) source for one code and two libraries.

1. [World Magnetic Anomaly Model (WMAM)](./apps/WMAM/README.md)
2. [Iterative Globel Model Fitting Library (globlibi)](./libs/globlibi/README.md)
3. [SLATEC Common Mathematical Library 4.1](./libs/slatec/README.md)

Example submission scripts can be found under [./tests/WMAM/submit](./tests/WMAM/submit).
Input data for the L=200, 300 and 720 cases (where L is the spherical harmonic degree)
can be found under [./tests/WMAM/shd](./tests/WMAM/shd).

The input data for L=1440 and L=2000 are too big to store within a git repo.

Please note, it is recommended that you first copy the `./tests/WMAM` folder to
`${HOME/home/work}/tests/WMAM`, submitting your jobs from a location outside the
eCSE04-8 repo directory.