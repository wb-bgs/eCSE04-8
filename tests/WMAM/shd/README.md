WMAM Input Data
===============

A WMAM run requires three input files, reference model (`coef_1990_15.dat.bin`),
observational data with coordinates (`wdmam_geocentric.dat.bin`) and starting model
(`model.in.bin`). The model files contain spherical harmonic coefficient values.

The input data differs according to the spherical harmonic degree and within this
folder you can find the input data sets for degrees [200](./200), [300](./300) and
[720](./720).

The starting model and observational data for degrees [1440](./1440) and [2000](./2000) are too big
to store within a git repo and so those folders just contain a `README.md` file that explains
how to generate the input data.

Generating input data requires running [various scripts](./scripts).
There is a script for creating the reference model file ([`coef_1990_15.dat`](./scripts/make_ref_model_file)), the observational data file ([`wdmam_geocentric.dat`](./scripts/make_wdmam_file)) and the starting model file ([`model.in`](./scripts/make_guess_model_file)). In addition, there is a script for converting these
three files from ASCII to binary ([`make_binary_files.sh`](./scripts/make_binary_files.sh)), allowing parallel processes to read the data. 