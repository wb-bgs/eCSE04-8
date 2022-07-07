#!/bin/bash --login

COMPILER_TAG=cce12
MPI_TAG=cmpich8
COMMS_TAG=ofi
SUBMIT_PRFX=${COMPILER_TAG}_${MPI_TAG}_${COMMS_TAG}

SUBMIT_PATH="./200/${COMPILER_TAG}/${COMMS_TAG}"
SBATCH_ARGS="--qos=lowpriority"

jobid=$(sbatch --parsable ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n01_c008.ll)

jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n01_c016.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n01_c032.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n01_c064.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n01_c128.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n02.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n04.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n08.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n16.ll)

sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n32.ll
