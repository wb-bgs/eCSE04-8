# ./plot_summary_data.sh 1440 0.1 1.6


ROOT=${HOME/home/work}
PYPP_HOME=${ROOT}/utils/pypp/3.9.4.1
SCRIPTS_HOME=${ROOT}/tests/scripts

APP_LABEL=WMAM
APP_VERSION=$3
APP_NAME="${APP_LABEL} ${APP_VERSION}"
DEGREE=$1
RESOLUTION=$2
TARGET_PATH=cce12/cmpich8-ofi
RESULTS_SUFFIX="${APP_VERSION}_scaling"
RESULTS_HOME=${ROOT}/tests/${APP_LABEL}/results.${RESULTS_SUFFIX}/${DEGREE}
PLOTS_HOME=${ROOT}/tests/${APP_LABEL}/plots.${RESULTS_SUFFIX}
SUMMARY_DATA=${RESULTS_HOME}/${TARGET_PATH}/summary.txt
SCRIPTS_HOME=${ROOT}/tests/${APP_LABEL}/scripts/plotting


TITLE="Strong Scaling ${APP_NAME} (ilg=${DEGREE}, deg_res=${RESOLUTION})"
DESC="ARCHER2 |Cray PE 21.09|Cray MPICH 8.1.9 (OFI)"

mkdir -p ${PLOTS_HOME}

. ${PYPP_HOME}/bin/activate

PLOT_NAME=${PLOTS_HOME}/wmam-${DEGREE}-runtime
python ${SCRIPTS_HOME}/plot_results.py ${SUMMARY_DATA} 0 "Core Count" 1 "Runtime (s)" 1.0 "${TITLE}" "${PLOT_NAME}" 0.5 "upper right" 0 0 "${DESC}" "0.12,0.94,-0.05"

PLOT_NAME=${PLOTS_HOME}/wmam-${DEGREE}-efficiency
python ${SCRIPTS_HOME}/plot_results.py ${SUMMARY_DATA} 0 "Core Count" 3 "Parallel Efficiency" 1.0 "${TITLE}" "${PLOT_NAME}" 0.5 "upper right" 0 0 "${DESC}" "0.6,0.94,-0.05"

. ${PYPP_HOME}/bin/deactivate
