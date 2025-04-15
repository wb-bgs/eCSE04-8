#!/bin/bash

# ./calculate_energy.sh gpu 10.2.0-nvfortran gnu10 ompi4 ucx 5.0 200 6619620

ROOT=${HOME/home/work}
PYPP_HOME=${ROOT}/utils/pypp/py38
SCRIPTS_HOME=${ROOT}/tests/WMAM/scripts

HARDWARE=$1
PE_RELEASE=$2
COMPILER_LABEL=$3
MPI_LABEL=$4
COMMS_LABEL=$5
APP_VERSION=$6
DEGREE=$7
JOBID=$8

APP_LABEL=WMAM
TARGET_PATH=${PE_RELEASE}/${COMPILER_LABEL}/${MPI_LABEL}-${COMMS_LABEL}
RESULTS_HOME=${ROOT}/tests/${APP_LABEL}/results/${DEGREE}
CALCULATE_SCRIPT=${SCRIPTS_HOME}/analysis/nvidia-smi-energy.py


set -

. ${PYPP_HOME}/bin/activate


python ${CALCULATE_SCRIPT} -p ${RESULTS_HOME}/${TARGET_PATH}/n\*/tpn\*/tpg\*/tpt\*/${JOBID}/nvsmi-power-\*.out


conda deactivate
