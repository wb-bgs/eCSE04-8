#!/bin/bash --login

COMPILER_TAG=cce12
MPI_TAG=cmpich8
COMMS_TAG=ofi
SUBMIT_PRFX=${COMPILER_TAG}_${MPI_TAG}_${COMMS_TAG}

SUBMIT_PATH="./2000/${COMPILER_TAG}/${COMMS_TAG}"
SBATCH_ARGS="--time=24:00:00 --qos=lowpriority"

jobid=$(sbatch --parsable ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n0128.ll)

jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n0256.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n0512.ll)

sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n1024.ll
