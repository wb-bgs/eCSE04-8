#!/bin/bash

# ./verify_output.sh 22.12 cce15 cmpich8 ofi 3.6 200 2784014
# ./verify_output.sh 10.2.0 gnu10 ompi4 ucx 5.0 200 2784014

ROOT=${HOME/home/work}
PYPP_HOME=${ROOT}/utils/pypp/py38
SCRIPTS_HOME=${ROOT}/tests/WMAM/scripts

PE_RELEASE=$1
COMPILER_LABEL=$2
MPI_LABEL=$3
COMMS_LABEL=$4
APP_VERSION=$5
DEGREE=$6
JOBID=$7

APP_LABEL=WMAM
TARGET_PATH=${PE_RELEASE}/${COMPILER_LABEL}/${MPI_LABEL}-${COMMS_LABEL}
RESULTS_HOME=${ROOT}/tests/${APP_LABEL}/results/${DEGREE}
VERIFY_SCRIPT=${SCRIPTS_HOME}/verification/verify_output.py

REFERENCE_HOME="${RESULTS_HOME}/ref/Results"

set -

. ${PYPP_HOME}/bin/activate

python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n*/tpn*/tpt*/${JOBID}/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out"

conda deactivate
