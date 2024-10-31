#!/bin/bash

#SBATCH --job-name=WMAM
#SBATCH --time=02:00:00
#SBATCH --exclusive
#SBATCH --nodes=256
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=16
#SBATCH --account=[budget code]
#SBATCH --partition=standard
#SBATCH --qos=standard


ulimit -c unlimited
ulimit -S -s unlimited


module -q load wmam/5.0

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

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OMP_PLACES=cores
export OMP_STACKSIZE=1G


srun ${SRUN_PARAMS} mod_wmam_020 ${WMAM_PARAMS}
