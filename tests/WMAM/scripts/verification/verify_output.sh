# ./verify_output.sh 200 888326


ROOT=${HOME/home/work}
PYPP_HOME=${ROOT}/utils/pypp/3.9.4.1
SCRIPTS_HOME=${ROOT}/tests/WMAM/scripts

APP_LABEL=WMAM
DEGREE=$1
TARGET_PATH=cce12/cmpich8-ofi
RESULTS_HOME=${ROOT}/tests/${APP_LABEL}/results/${DEGREE}
VERIFY_SCRIPT=${SCRIPTS_HOME}/verification/verify_output.py

REFERENCE_DEGREE=""
if [ ${DEGREE} -eq 200 ]; then
  REFERENCE_HOME="~/../shared/arc/eCSE_medium/Results"
elif [ ${DEGREE} -eq 300 ]; then
  REFERENCE_HOME="~/../shared/arc/eCSE_long/Results"
elif [ ${DEGREE} -eq 720 ]; then
  REFERENCE_HOME="~/../shared/arc/720/Results"
  VERIFY_SCRIPT=${SCRIPTS_HOME}/verification/verify_output2.py
elif [ ${DEGREE} -eq 1440 ]; then
  REFERENCE_HOME="~/../shared/arc/mod_wdmam_4Nick_1024/Results"
elif [ ${DEGREE} -eq 2000 ]; then
  REFERENCE_HOME="~/../shared/arc/eCSE_large_problem/eCSE2021-v1.3/Results"
else
  echo "Error, unrecognised degree value."
  exit 0
fi

. ${PYPP_HOME}/bin/activate

python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n1/c64/$2/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out"

#python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n1/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n1.out
#python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n2/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n2.out
#python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n4/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n4.out
#python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n8/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n8.out
#python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n16/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n16.out
#python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n32/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n32.out
#python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n64/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n64.out
#python ${VERIFY_SCRIPT} "${RESULTS_HOME}/${TARGET_PATH}/n128/c*/*/Results/model_No_P.out" "${REFERENCE_HOME}/model_No_P.out" &> comp_${DEGREE}_n128.out

. ${PYPP_HOME}/bin/deactivate
