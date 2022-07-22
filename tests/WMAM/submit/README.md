WMAM Job Submission
===================

This folder contains contains numerous Slurm submission scripts
for running specific versions of WMAM over specific numbers of nodes
for specific spherical harmonic degrees.

There are also submission scripts for running WMAM with various profiling
and debugging tools, see the [tests folder](./tests).

The WMAM application itself is called with four arguments, the spherical harmonic degree (L),
the resolution degree (d), the inversion scheme id (Polak-Ribiere=1, Conjugate-Gradient=2)
and the damping factor (usually 5.0).

The spherical harmonic degree and resolution degree values are paired:
the following combinations (L,d) have been tested so far, (200,1.0), (300,0.5),
(720,0.25), (1440,0.1) and (2000,0.05).

Finally, there is a set of launch scripts (for L=200, 300, 720, 1440 and 2000)
that submit the runs required to measure strong scaling performance. In general,
scaling runs will be submitted with `--qos=lowpriority` to save on resource use.