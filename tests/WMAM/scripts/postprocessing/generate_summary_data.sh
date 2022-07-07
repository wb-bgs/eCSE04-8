# ./generate_summary_data.sh 1440 1.6


ROOT=${HOME/home/work}
APP_LABEL=WMAM
APP_VERSION=$2
DEGREE=$1
RESULTS_SUFFIX="${APP_VERSION}_scaling"
TARGET_PATH=cce12/cmpich8-ofi
PYPP_HOME=${ROOT}/utils/pypp/3.9.4.1
SCRIPTS_HOME=${ROOT}/tests/${APP_LABEL}/scripts/postprocessing
RESULTS_HOME=${ROOT}/tests/${APP_LABEL}/results.${RESULTS_SUFFIX}/${DEGREE}


. ${PYPP_HOME}/bin/activate

python ${SCRIPTS_HOME}/get_performance_metric.py ${APP_LABEL} "${RESULTS_HOME}/${TARGET_PATH}/n*/c*/*/${APP_LABEL}.o" > ${RESULTS_HOME}/${TARGET_PATH}/summary.txt

. ${PYPP_HOME}/bin/deactivate
