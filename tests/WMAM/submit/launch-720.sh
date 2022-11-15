#!/bin/bash

SUBMIT_PATH="./720"
SBATCH_ARGS="--time=12:00:00 --qos=lowpriority"

jobid=$(sbatch --parsable ${SBATCH_ARGS} ${SUBMIT_PATH}/n002.ll)

jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n004.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n008.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n016.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n032.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n064.ll)

sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n128.ll
