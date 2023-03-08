#!/bin/bash

# make_ascii_files.sh 200 21.09 cce12 cmpich8-ofi 1 64 3209856

SHD=$1
CPE_VERSION=$2
COMPILER_LABEL=$3
MPI_LABEL=$4
NNODES=$5
NCORES=$6
JOBID=$7

SCRIPTS_PATH=${HOME/home/work}/tests/WMAM/shd/scripts
RESULTS_PATH=${HOME/home/work}/tests/WMAM/results/${SHD}/${CPE_VERSION}/${COMPILER_LABEL}/${MPI_LABEL}/n${NNODES}/c${NCORES}/${JOBID}/Results

module -q load cray-python

python ${SCRIPTS_PATH}/make_ascii_file.py ${RESULTS_PATH} fit_No_P.out.bin
python ${SCRIPTS_PATH}/make_ascii_file.py ${RESULTS_PATH} fit_damp.out.bin