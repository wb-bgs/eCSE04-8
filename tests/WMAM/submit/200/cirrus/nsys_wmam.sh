#!/bin/bash

MPI_RANKS_PER_NODE=$1
APP_EXE=$2
APP_PARAMS=$3

NSYS_PARAMS="profile --trace=mpi,openmp,cuda --mpi-impl=openmpi --cuda-memory-usage=true --output=wmam-baseline-${SLURM_PROCID} --env-var=NSYS_MPI_STORE_TEAMS_PER_RANK=1"

nsys ${NSYS_PARAMS} ${APP_EXE} ${APP_PARAMS}
