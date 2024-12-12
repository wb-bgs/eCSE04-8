#!/bin/bash

# ./verify_output.sh cpu 10.2.0 gnu10 ompi4 ucx 5.0 200 6619617
# ./verify_output.sh gpu 10.2.0-nvfortran gnu10 ompi4 ucx 5.0 200 6619620

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
VERIFY_SCRIPT=${SCRIPTS_HOME}/verification/verify_output.py

REFERENCE_HOME="${RESULTS_HOME}/ref/Results"

set -

. ${PYPP_HOME}/bin/activate

if [[ "$HARDWARE" == "cpu" ]]; then
  python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n*/tpn*/tpt*/${JOBID}/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out"
fi

if [[ "$HARDWARE" == "gpu" ]]; then
  python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n*/tpn*/tpg*/tpt*/${JOBID}/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out"
fi

conda deactivate
