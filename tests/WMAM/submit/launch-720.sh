#!/bin/bash --login

COMPILER_TAG=cce12
MPI_TAG=cmpich8
COMMS_TAG=ofi
SUBMIT_PRFX=${COMPILER_TAG}_${MPI_TAG}_${COMMS_TAG}

SUBMIT_PATH="./720/${COMPILER_TAG}/${COMMS_TAG}"
SBATCH_ARGS="--time=12:00:00 --qos=lowpriority"

jobid=$(sbatch --parsable ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n002.ll)

jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n004.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n008.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n016.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n032.ll)
jobid=$(sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n064.ll)

sbatch --parsable --dependency=afterok:${jobid} ${SBATCH_ARGS} ${SUBMIT_PATH}/${SUBMIT_PRFX}_n128.ll
