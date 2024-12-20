#!/bin/bash

#SBATCH --job-name=WMAM
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --nodes=512
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --account=ecsead08
#SBATCH --partition=standard
#SBATCH --qos=standard


ulimit -c unlimited
ulimit -S -s unlimited


module -q restore
module -q load PrgEnv-cray


# setup resource-related environment
NNODES=${SLURM_JOB_NUM_NODES}
NTASKSPN=`echo "${SLURM_TASKS_PER_NODE}" | cut -d'(' -f1`
NCORESPT=${SLURM_CPUS_PER_TASK}

NCORESPN=`expr ${NTASKSPN} \* ${NCORESPT}`
NCORES=`expr ${NNODES} \* ${NCORESPN}`
NTASKS=`expr ${NNODES} \* ${NTASKSPN}`


PE_NAME=${PE_MPICH_FIXED_PRGENV}
PE_RELEASE=22.12

ROOT=${HOME/home/work}
DEGREE=1440
RESOLUTION=0.1
SCHEME=1
DAMPFAC=5.0
SERIALRD=1
APP_NAME=WMAM
APP_VERSION=4.1
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

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_STACKSIZE=1G


# setup app run directory and input folder
mkdir -p ${APP_RUN_PATH}/Data

# copy input data
cp ${APP_RUN_ROOT}/shd/${DEGREE}/* ${APP_RUN_PATH}/Data/

# setup output folder
mkdir ${APP_RUN_PATH}/Results
lfs setstripe -c -1 ${APP_RUN_PATH}/Results


# run app
RUN_START=$(date +%s.%N)
echo -e "Launching ${APP_EXE_NAME} ${APP_VERSION} (${PE_RELEASE},${APP_COMPILER_LABEL},${APP_MPI_LABEL}-${APP_COMMS_LABEL}) with l=${DEGREE} across ${NNODES} node(s), ${NTASKSPN} task(s) per node and ${NCORESPT} thread(s) per task.\n" > ${APP_OUTPUT}

srun ${SRUN_PARAMS} ${APP_EXE} ${APP_PARAMS} &>> ${APP_OUTPUT}

RUN_STOP=$(date +%s.%N)
RUN_TIME=$(echo "${RUN_STOP} - ${RUN_START}" | bc)
echo -e "\nsrun time: ${RUN_TIME}" >> ${APP_OUTPUT}


# tidy up
mkdir ${APP_RUN_PATH}/${SLURM_JOB_ID}
mv ${APP_OUTPUT} ${APP_RUN_PATH}/${SLURM_JOB_ID}/
mv ${APP_RUN_PATH}/Results ${APP_RUN_PATH}/${SLURM_JOB_ID}/
rm -rf ${APP_RUN_PATH}/Data
