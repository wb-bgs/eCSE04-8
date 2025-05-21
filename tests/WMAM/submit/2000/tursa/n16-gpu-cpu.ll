#!/bin/bash

#SBATCH --job-name=WMAM
#SBATCH --time=00:20:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --account=z01
#SBATCH --partition=gpu-a100-40
#SBATCH --qos=dev
#SBATCH --gres=gpu:4


cat $0


ulimit -c unlimited
ulimit -S -s unlimited


module load gcc/12.2.0
module load openmpi/4.1.5-gcc12-cpu


# setup resource-related environment
NNODES=${SLURM_JOB_NUM_NODES}
NCORESPN=${SLURM_CPUS_ON_NODE}
NTASKSPN=4
NCORESPT=8

NTASKS=`expr ${NNODES} \* ${NTASKSPN}`
NCORES=`expr ${NNODES} \* ${NTASKSPN} \* ${NCORESPT}`


PE_NAME=GNU
PE_RELEASE=12.2.0

ROOT=${HOME}
DEGREE=2000
RESOLUTION=0.05
SCHEME=1
DAMPFAC=5.0
SERIALRD=0
APP_NAME=WMAM
APP_VERSION=5.0
APP_EXE_NAME=mod_wmam_020
APP_EXE=${ROOT}/apps/${APP_NAME}/${APP_VERSION}/${PE_NAME}/${PE_RELEASE}/release/bin/${APP_EXE_NAME}
APP_MPI_LABEL=ompi4
APP_COMMS_LABEL=ucx
APP_COMPILER_LABEL=gnu10
APP_RUN_ROOT=${ROOT}/tests/${APP_NAME}
APP_RUN_PATH=${APP_RUN_ROOT}/results/${DEGREE}/${PE_RELEASE}/${APP_COMPILER_LABEL}/${APP_MPI_LABEL}-${APP_COMMS_LABEL}/n${NNODES}/tpn${NTASKSPN}/tpt${NCORESPT}/${SLURM_JOB_ID}
APP_PARAMS="${DEGREE} ${RESOLUTION} ${SCHEME} ${DAMPFAC} ${SERIALRD}"


SRUN_PARAMS="--nodes=${NNODES} --ntasks-per-node=${NTASKSPN} --cpus-per-task=${NCORESPT} --hint=nomultithread --distribution=block:block --unbuffered --chdir=${APP_RUN_PATH}"

export OMP_NUM_THREADS=${NCORESPT}
export OMP_PLACES=cores


# setup app run directory and input folder
mkdir -p ${APP_RUN_PATH}/Data

# copy input data
cp ${APP_RUN_ROOT}/shd/${DEGREE}/* ${APP_RUN_PATH}/Data/

# setup output folder
mkdir ${APP_RUN_PATH}/Results
lfs setstripe -c -1 ${APP_RUN_PATH}/Results


# run app
RUN_START=$(date +%s.%N)
echo -e "\n\nLaunching ${APP_EXE_NAME} ${APP_VERSION} (${PE_RELEASE},${APP_COMPILER_LABEL},${APP_MPI_LABEL}-${APP_COMMS_LABEL}) with l=${DEGREE} across ${NNODES} node(s), ${NTASKSPN} task(s) per node and ${OMP_NUM_THREADS} thread(s) per task.\n"

srun ${SRUN_PARAMS} ${APP_EXE} ${APP_PARAMS}

RUN_STOP=$(date +%s.%N)
RUN_TIME=$(echo "${RUN_STOP} - ${RUN_START}" | bc)
echo -e "\nsrun time: ${RUN_TIME}"


# tidy up
mv ./slurm-${SLURM_JOB_ID}.out ${APP_RUN_PATH}/${APP_NAME}.o
rm -rf ${APP_RUN_PATH}/Data
