#!/bin/bash

#SBATCH --job-name=WMAM
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --account=ecsead08
#SBATCH --partition=standard
#SBATCH --qos=standard


ulimit -c unlimited


module -q restore
module -q load PrgEnv-cray


# setup resource-related environment
NNODES=${SLURM_JOB_NUM_NODES}
NTASKSPN=`echo "${SLURM_TASKS_PER_NODE}" | cut -d'(' -f1`
NCORES=`expr ${NTASKSPN} \* ${SLURM_CPUS_PER_TASK}`
if [ "${NNODES}" -gt "1" ]; then
  NCORES=`expr ${NNODES} \* ${NCORES}`
fi
export OMP_NUM_THREADS=1


PE_NAME=${PE_MPICH_FIXED_PRGENV}
PE_RELEASE=22.12

ROOT=${HOME/home/work}
DEGREE=1440
RESOLUTION=0.1
SCHEME=1
DAMPFAC=5.0
SERIALRD=0
APP_NAME=WMAM
APP_VERSION=4.0
APP_EXE_NAME=mod_wmam_020
APP_EXE=${ROOT}/apps/${APP_NAME}/${APP_VERSION}/${PE_NAME}/${PE_RELEASE}/release/bin/${APP_EXE_NAME}
APP_MPI_LABEL=cmpich8
APP_COMMS_LABEL=ofi
APP_COMPILER_LABEL=cce15
APP_RUN_ROOT=${ROOT}/tests/${APP_NAME}
APP_RUN_PATH=${APP_RUN_ROOT}/results/${DEGREE}/${PE_RELEASE}/${APP_COMPILER_LABEL}/${APP_MPI_LABEL}-${APP_COMMS_LABEL}/n${NNODES}/c${NCORES}
APP_PARAMS="${DEGREE} ${RESOLUTION} ${SCHEME} ${DAMPFAC} ${SERIALRD}"
APP_OUTPUT=${APP_RUN_PATH}/${APP_NAME}.o

export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
SRUN_PARAMS="--distribution=block:block --hint=nomultithread --unbuffered --chdir=${APP_RUN_PATH}"


# setup app run directory and input folder
mkdir -p ${APP_RUN_PATH}/Data

# copy input data
cp ${APP_RUN_ROOT}/shd/${DEGREE}/* ${APP_RUN_PATH}/Data/

# setup output folder
mkdir ${APP_RUN_PATH}/Results
lfs setstripe -c -1 ${APP_RUN_PATH}/Results


# run app
RUN_START=$(date +%s.%N)
echo -e "Launching ${APP_EXE_NAME} ${APP_VERSION} (${PE_RELEASE},${APP_COMPILER_LABEL},${APP_MPI_LABEL}-${APP_COMMS_LABEL}) with L=${DEGREE} over ${NCORES} core(s) across ${NNODES} node(s).\n" > ${APP_OUTPUT}

srun ${SRUN_PARAMS} ${APP_EXE} ${APP_PARAMS} &>> ${APP_OUTPUT}

RUN_STOP=$(date +%s.%N)
RUN_TIME=$(echo "${RUN_STOP} - ${RUN_START}" | bc)
echo -e "\nsrun time: ${RUN_TIME}" >> ${APP_OUTPUT}


# tidy up
mkdir ${APP_RUN_PATH}/${SLURM_JOB_ID}
mv ${APP_OUTPUT} ${APP_RUN_PATH}/${SLURM_JOB_ID}/
mv ${APP_RUN_PATH}/Results ${APP_RUN_PATH}/${SLURM_JOB_ID}/
rm -rf ${APP_RUN_PATH}/Data
