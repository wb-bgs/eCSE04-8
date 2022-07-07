#!/bin/bash

SHD=$1

module -q load cray-python

python make_binary_file.py ${HOME/home/work}/tests/WMAM/shd/${SHD} coef_1990_15.dat
python make_binary_file.py ${HOME/home/work}/tests/WMAM/shd/${SHD} model.in
python make_binary_file.py ${HOME/home/work}/tests/WMAM/shd/${SHD} wdmam_geocentric.dat
