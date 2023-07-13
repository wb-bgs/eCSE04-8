#!/bin/bash

# ./verify_output.sh 22.12 cce15 3.6 200 2784014

ROOT=${HOME/home/work}
PYPP_HOME=${ROOT}/utils/pypp
SCRIPTS_HOME=${ROOT}/tests/WMAM/scripts

CPE_RELEASE=$1
COMPILER_LABEL=$2
APP_VERSION=$3
DEGREE=$4
JOBID=$5

APP_LABEL=WMAM
TARGET_PATH=${CPE_RELEASE}/${COMPILER_LABEL}/cmpich8-ofi
RESULTS_HOME=${ROOT}/tests/${APP_LABEL}/results/${DEGREE}
VERIFY_SCRIPT=${SCRIPTS_HOME}/verification/verify_output.py

REFERENCE_DEGREE=""
if [ ${DEGREE} -eq 200 ]; then
  REFERENCE_HOME="~/../shared/arc/eCSE_medium/Results"
elif [ ${DEGREE} -eq 300 ]; then
  REFERENCE_HOME="~/../shared/arc/eCSE_long/Results"
elif [ ${DEGREE} -eq 720 ]; then
# REFERENCE_HOME="~/../shared/arc/720/Results"
  REFERENCE_HOME="${HOME/home/work}/tests/WMAM/results/720/22.12/cce15/cmpich8-ofi/n2/c256/2528530/Results"
elif [ ${DEGREE} -eq 1440 ]; then
# REFERENCE_HOME="~/../shared/arc/mod_wdmam_4Nick_1024/Results"
  REFERENCE_HOME="${HOME/home/work}/tests/WMAM/results/scaling/3.0/1440/cce15/cmpich8-ofi/n128/c16384/2034985/Results"
elif [ ${DEGREE} -eq 2000 ]; then
# REFERENCE_HOME="~/../shared/arc/eCSE_large_problem/eCSE2021-v1.3/Results"
  REFERENCE_HOME="${HOME/home/work}/tests/WMAM/results/scaling/3.1/2000/cce15/cmpich8-ofi/n256/c32768/2052655/Results"
else
  echo "Error, unrecognised degree value."
  exit 0
fi

. ${PYPP_HOME}/bin/activate

python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n*/c*/${JOBID}/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out"

deactivate
