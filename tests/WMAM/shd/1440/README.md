WMAM input data for *L* = 1440
==============================

1. Generate reference model

```bash
# requires "bggm2017.dat"
# outputs "coef_1990_15.dat"

../scripts/make_ref_model_file
```

2. Generate guess model

```bash
# requires "model_No_P.out"
# outputs "model.in" (and "model.tmp") 

../scripts/make_guess_model_file 1440
```

3. Generate input data

```bash
# requires "wdmam.asc" and two utility codes, "gmt" and "gd2gc_cmdln" (see the top-level utils folder)
# outputs "wdmam_geocentric.dat" (and "wdmam0p2_gdlo-gdla-dF.tmp",
#                                     "wdmam0p2_gdlo-gdla-id.tmp",
#                                     "wdmam_geocentric.tmp",
#                                     "wdmam_geodetic.tmp") 

../scripts/make_wdmam_file 0.1
```

The `wdmam_geocentric.dat` and `model.in` files are too big to store in this repository.

Those two files can of course be generated by running the `make_wdmam_file` and `make_guess_model_file` scripts.
However, those scripts require precursor data files that are also too large to store here.
Therfore, instructions for downloading the precursor data files can be found at [../data/README.md](../data/README.md).
