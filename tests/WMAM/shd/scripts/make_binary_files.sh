#!/bin/bash

SHD=$1
SHD_PATH=${HOME/home/work}/tests/WMAM/shd/${SHD}

module -q load cray-python

python make_binary_file.py ${SHD_PATH} coef_1990_15.dat
python make_binary_file.py ${SHD_PATH} model.in
python make_binary_file.py ${SHD_PATH} wdmam_geocentric.dat
