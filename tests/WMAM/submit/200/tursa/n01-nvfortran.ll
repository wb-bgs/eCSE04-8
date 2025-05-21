#!/bin/bash

#SBATCH --job-name=WMAM
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --account=z01
#SBATCH --partition=gpu-a100-40
#SBATCH --qos=dev
#SBATCH --gres=gpu:4
#SBATCH --gpu-freq=1000


cat $0


ulimit -c unlimited
ulimit -S -s unlimited


module load gcc/12.2.0
module load nvhpc/23.5-nompi
module load openmpi/4.1.5-nvhpc235-cuda12


# setup resource-related environment
NNODES=${SLURM_JOB_NUM_NODES}
NCORESPN=`expr ${SLURM_CPUS_ON_NODE} \/ 2`
NGPUSPN=${SLURM_GPUS_ON_NODE}

NTASKSPN=${NGPUSPN}
NTASKSPG=`expr ${NTASKSPN} \/ ${NGPUSPN}`
NCORESPT=`expr ${NCORESPN} \/ ${NTASKSPN}`

NCORES=`expr ${NNODES} \* ${NCORESPN}`
NTASKS=`expr ${NNODES} \* ${NTASKSPN}`

export OMP_NUM_THREADS=1
export OMP_PLACES=cores


PE_NAME=GNU
PE_RELEASE=12.2.0-nvfortran

ROOT=${HOME}
DEGREE=200
RESOLUTION=1.0
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
APP_RUN_PATH=${APP_RUN_ROOT}/results/${DEGREE}/${PE_RELEASE}/${APP_COMPILER_LABEL}/${APP_MPI_LABEL}-${APP_COMMS_LABEL}/n${NNODES}/tpn${NTASKSPN}/tpg${NTASKSPG}/tpt${OMP_NUM_THREADS}/${SLURM_JOB_ID}
APP_PARAMS="${DEGREE} ${RESOLUTION} ${SCHEME} ${DAMPFAC} ${SERIALRD}"


SRUN_PARAMS="--nodes=${NNODES} --tasks-per-node=${NTASKSPN} --cpus-per-task=${NCORESPT} --hint=nomultithread --distribution=block:block --unbuffered --chdir=${APP_RUN_PATH}"


# setup app run directory and input folder
mkdir -p ${APP_RUN_PATH}/Data

# copy input data
cp ${APP_RUN_ROOT}/shd/${DEGREE}/* ${APP_RUN_PATH}/Data/

# setup output folder
mkdir ${APP_RUN_PATH}/Results
lfs setstripe -c -1 ${APP_RUN_PATH}/Results


# run app
RUN_START=$(date +%s.%N)
echo -e "\n\nLaunching ${APP_EXE_NAME} ${APP_VERSION} (${PE_RELEASE},${APP_COMPILER_LABEL},${APP_MPI_LABEL}-${APP_COMMS_LABEL}) with l=${DEGREE} across ${NNODES} node(s), ${NTASKSPN} task(s) per node, ${NTASKSPG} task(s) per GPU and ${OMP_NUM_THREADS} thread(s) per task.\n"

srun ${SRUN_PARAMS} \
    ${APP_RUN_ROOT}/submit/${DEGREE}/tursa/wmam.sh ${NTASKSPN} 1 "${APP_RUN_PATH}" "${APP_EXE}" "${APP_PARAMS}"

RUN_STOP=$(date +%s.%N)
RUN_TIME=$(echo "${RUN_STOP} - ${RUN_START}" | bc)
echo -e "\nsrun time: ${RUN_TIME}"


# tidy up
mv ./slurm-${SLURM_JOB_ID}.out ${APP_RUN_PATH}/${APP_NAME}.o
rm -rf ${APP_RUN_PATH}/Data
