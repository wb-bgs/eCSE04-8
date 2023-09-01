#!/bin/bash

#SBATCH --job-name=WMAM
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --account=[budget code]
#SBATCH --partition=standard
#SBATCH --qos=standard


ulimit -c unlimited


module -q load wmam/4.0

DEGREE=200
RESOLUTION=1.0
SCHEME=1
DAMPFAC=5.0
SERIALRD=0

WMAM_PARAMS="${DEGREE} ${RESOLUTION} ${SCHEME} ${DAMPFAC} ${SERIALRD}"

mkdir ${SLURM_SUBMIT_DIR}/Results
lfs setstripe -c -1 ${SLURM_SUBMIT_DIR}/Results


export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}
SRUN_PARAMS="--distribution=block:block --hint=nomultithread --unbuffered"


srun ${SRUN_PARAMS} mod_wmam_020 ${WMAM_PARAMS}
