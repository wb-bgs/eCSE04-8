#!/bin/bash

SUBMIT_PATH="./200"
SBATCH_ARGS="--qos=lowpriority"

jobid=$(sbatch --parsable ${SBATCH_ARGS} ${SUBMIT_PATH}/n01_c008.ll)

jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n01_c016.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n01_c032.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n01_c064.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n01.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n02.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n04.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n08.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n16.ll)

sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/n32.ll
