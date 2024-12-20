#!/bin/bash

#SBATCH --job-name=WMAM
#SBATCH --time=00:20:00
#SBATCH --exclusive
#SBATCH --nodes=2
#SBATCH --account=z04
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --gres=gpu:4


cat $0


ulimit -c unlimited
ulimit -S -s unlimited


module -s load nvidia/nvhpc-nompi/24.5
module -s load openmpi/4.1.6-cuda-12.4-nvfortran


# setup resource-related environment
NNODES=${SLURM_JOB_NUM_NODES}
NCORESPN=${SLURM_CPUS_ON_NODE}
NGPUSPN=${SLURM_GPUS_ON_NODE}

NTASKSPN=${NGPUSPN}
NTASKSPG=`expr ${NTASKSPN} \/ ${NGPUSPN}`
NCORESPT=`expr ${NCORESPN} \/ ${NTASKSPN}`

NCORES=`expr ${NNODES} \* ${NCORESPN}`
NTASKS=`expr ${NNODES} \* ${NTASKSPN}`

export OMP_NUM_THREADS=1
export OMP_PLACES=cores


PE_NAME=GNU
PE_RELEASE=10.2.0-nvfortran

ROOT=${HOME/home/work}
DEGREE=300
RESOLUTION=0.5
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
APP_RUN_PATH=${APP_RUN_ROOT}/results/${DEGREE}/${PE_RELEASE}/${APP_COMPILER_LABEL}/${APP_MPI_LABEL}-${APP_COMMS_LABEL}/n${NNODES}/tpn${NTASKSPN}/tpg${NTASKSPG}/tpt${OMP_NUM_THREADS}
APP_PARAMS="${DEGREE} ${RESOLUTION} ${SCHEME} ${DAMPFAC} ${SERIALRD}"


SRUN_PARAMS="--ntasks=${NTASKS} --tasks-per-node=${NTASKSPN} --cpus-per-task=${NCORESPT} --hint=nomultithread --unbuffered --chdir=${APP_RUN_PATH}"


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

export HOME=${HOME/home/work}


srun ${SRUN_PARAMS} \
    ${APP_RUN_ROOT}/submit/${DEGREE}/cirrus/nsys_wmam.sh ${NTASKSPN} "${APP_EXE}" "${APP_PARAMS}"


RUN_STOP=$(date +%s.%N)
RUN_TIME=$(echo "${RUN_STOP} - ${RUN_START}" | bc)
echo -e "\nsrun time: ${RUN_TIME}"


# tidy up
mkdir ${APP_RUN_PATH}/${SLURM_JOB_ID}
mv ./slurm-${SLURM_JOB_ID}.out ${APP_RUN_PATH}/${SLURM_JOB_ID}/${APP_NAME}.o
mv ${APP_RUN_PATH}/Results ${APP_RUN_PATH}/${SLURM_JOB_ID}/
mv ${APP_RUN_PATH}/wmam-baseline-*.nsys-rep ${APP_RUN_PATH}/${SLURM_JOB_ID}/
rm -rf ${APP_RUN_PATH}/Data
