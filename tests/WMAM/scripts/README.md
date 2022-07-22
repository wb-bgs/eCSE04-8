WMAM Output Data
================

This folder contains contains various scripts for processing WMAM output data.

There is a [verification script](./verification/verify_output.sh) that reports on the differences between
a `model_No_P.out` output file and a suitable reference `model_No_P.out`.

There is a script for gathering performance data from scaling runs, see [`generate_summary_data.sh`](./postprocessing/generate_summary_data.sh)
and there is a script for plotting performance data ([`plot_summary_data.sh`](./plotting/plot_summary_data.sh)).

There is also a python script, [`plot_spectra.py`](./analysis/plot_spectra.py), for plotting power spectra.

Lastly, there is a rudimentary script for summing up memory use during program execution ([`memory_requirement.py`](./misc/memory_requirement.py))
and one other miscellaenous script that plots, amongst other things, the relationship between spherical harmonic degree and data point counts
([`plot_misc.py`](./plotting/plot_misc.py)).