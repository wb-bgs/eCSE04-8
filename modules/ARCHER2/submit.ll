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


module -q load wmam/3.7

DEGREE=200
RESOLUTION=1.0
SCHEME=1
DAMPFAC=5.0

mkdir ${SLURM_SUBMIT_DIR}/Results
lfs setstripe -c -1 ${SLURM_SUBMIT_DIR}/Results

srun --distribution=block:block --hint=nomultithread \
    mod_wmam_020 ${DEGREE} ${RESOLUTION} ${SCHEME} ${DAMPFAC}
