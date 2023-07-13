WMAM Precusor Data
==================

This readme file explains how to obtain the precursor data files required to produce the WMAM input files.

These precusor files should be placed in this `data` folder as that's the path referenced within the scripts
that generate the input files, see below.

1. [make_ref_model_file](../scripts/make_ref_model_file) takes `bggm2017.dat` and outputs `coef_1990_15.dat`.
2. [make_wdmam_file](../scripts/make_wdmam_file) takes `wdmam.asc` and outputs `wdmam_geocentric.dat.bin`
3. [make_guess_model_file](../scripts/make_guess_model_file) takes `model_No_P.out` and outputs `model.in.bin`.

The [bggm2017.dat](./bggm2017.dat) precusor file can be found in this folder. The other two precusor files
exceed GitHub's 100 MB file limit and so have to be obtained from outside this repo.

The `wdmam.asc` file can be downloaded from [http://wdmam.org/file/wdmam.asc](http://wdmam.org/file/wdmam.asc).

The `model_No_P.out` file can be downloaded using the University of Edinburgh DataSync service, click
[https://datasync.ed.ac.uk/index.php/s/kdpD1POVH5iMwBp](https://datasync.ed.ac.uk/index.php/s/kdpD1POVH5iMwBp)
- the password is "ecse04-8-precursor-data".
