WMAM input data for L=300
=========================

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

../scripts/make_guess_model_file 300
```

3. Generate input data

```bash
# requires "wdmam.asc", "gmt" and "gd2gc_cmdln"
# outputs "wdmam_geocentric.dat" (and "wdmam0p2_gdlo-gdla-dF.tmp",
#                                     "wdmam0p2_gdlo-gdla-id.tmp",
#                                     "wdmam_geocentric.tmp",
#                                     "wdmam_geodetic.tmp") 

../scripts/make_wdmam_file 0.5
```
