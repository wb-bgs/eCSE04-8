WMAM Job Submission
===================

This folder contains contains numerous Slurm submission scripts
for running specific versions of WMAM over specific numbers of nodes
for specific spherical harmonic degrees.

There are also submission scripts for running WMAM with various profiling
and debugging tools, see the [tests folder](./tests).

The WMAM code itself is called with four arguments, the spherical harmonic degree (*L*),
the resolution degree (*R*), the inversion scheme id (Polak-Ribiere=1, Conjugate-Gradient=2)
and the damping factor (usually 5.0).

The spherical harmonic degree and resolution degree values are paired:
the following combinations (*L*, *R*) have been tested so far, (200,1.0), (300,0.5),
(720,0.25), (1440,0.1) and (2000,0.05).


The input data files reside in a folder called `Data` and the output files
are written to folder called `Results`.

Three input files are required by WMAM.

1. `coef_1990_15.dat.bin`
The model coefficients consistent with the 1990 terrestrial core field.

2. `wdmam_geocentric.dat.bin`
The magnetic anomaly data.

3. `model.in.bin`
The starting coefficients for the output field mode.

Four files are produced by WMAM.

1. `WMAM.o`
Standard ouput/error data.

2. `model_No_P.out`
The final field model coefficients.

3. `fit_No_P.out.bin`
Diagnostic info relating to the geographic points covered by the anomaly data.

4. `fit_damp.out.bin`
Diagnostic info relating to the nodal points of the spherical harmonic model.

Those files that have a `.bin` suffix are read/written by MPI File I/O subroutines.
There are scripts for converting [ASCII files to binary input files](https://github.com/wb-bgs/eCSE04-8/blob/main/tests/WMAM/shd/scripts/make_binary_files.sh) and
for [converting binary output files to ASCII files](https://github.com/wb-bgs/eCSE04-8/blob/main/tests/WMAM/shd/scripts/make_ascii_files.sh).


Finally, there is a set of launch scripts (for *L* = 200, 300, 720, 1440 and 2000)
that submit the runs required to measure strong scaling performance. In general,
scaling runs will be submitted with `--qos=lowpriority` to save on resource use.
