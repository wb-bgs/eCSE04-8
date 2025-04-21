#!/bin/bash

# ./verify_scaling_output.sh 22.12 cce15 3.6 200

ROOT=${HOME/home/work}
PYPP_HOME=${ROOT}/utils/pypp/py38
SCRIPTS_HOME=${ROOT}/tests/WMAM/scripts

CPE_RELEASE=$1
COMPILER_LABEL=$2
APP_VERSION=$3
DEGREE=$4

APP_LABEL=WMAM
TARGET_PATH=${CPE_RELEASE}/${COMPILER_LABEL}/cmpich8-ofi
RESULTS_HOME=${ROOT}/tests/${APP_LABEL}/results/scaling/${APP_VERSION}/${DEGREE}
VERIFY_SCRIPT=${SCRIPTS_HOME}/verification/verify_output.py

REFERENCE_DEGREE=""
if [ ${DEGREE} -eq 200 ]; then
  REFERENCE_HOME="~/../shared/arc/eCSE_medium/Results"
elif [ ${DEGREE} -eq 300 ]; then
  REFERENCE_HOME="~/../shared/arc/eCSE_long/Results"
elif [ ${DEGREE} -eq 720 ]; then
# REFERENCE_HOME="~/../shared/arc/720/Results"
  REFERENCE_HOME="${HOME/home/work}/tests/WMAM/results//720/cce15/cmpich8-ofi/n2/c256/2528530/Results"
elif [ ${DEGREE} -eq 1440 ]; then
# REFERENCE_HOME="~/../shared/arc/mod_wdmam_4Nick_1024/Results"
  REFERENCE_HOME="${HOME/home/work}/tests/WMAM/results/scaling/3.0/1440/cce15/cmpich8-ofi/n128/c16384/2034985/Results"
elif [ ${DEGREE} -eq 2000 ]; then
# REFERENCE_HOME="~/../shared/arc/eCSE_large_problem/eCSE2021-v1.3/Results"
  REFERENCE_HOME="${HOME/home/work}/tests/WMAM/results.test/2000/cce15/cmpich8-ucx/n1024/c32768/1438712/Results"
else
  echo "Error, unrecognised degree value."
  exit 0
fi

. ${PYPP_HOME}/bin/activate

declare -a node_count=( "1" "2" "4" "8" "16" "32" "64" "128" )
for nc in "${node_count[@]}"; do
  python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n${nc}/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n${nc}.out
done

deactivate
