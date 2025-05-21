#!/bin/bash

MPI_RANKS_PER_NODE=$1
NVSMI_FREQUENCY=$2
APP_RUN_PATH=$3
APP_EXE=$4
APP_PARAMS=$5


# Bind the process to the correct GPU and NIC
export CUDA_VISIBLE_DEVICES=${SLURM_LOCALID}
export UCX_NET_DEVICES=mlx5_${SLURM_LOCALID}:1

if ! ((${SLURM_LOCALID} % ${MPI_RANKS_PER_NODE})); then

  # one MPI task per node will gather GPU power readings
  nvidia-smi --query-gpu=index,timestamp,power.draw --format=csv --loop=${NVSMI_FREQUENCY} &> ${APP_RUN_PATH}/nvsmi-power-${SLURM_PROCID}.out &

  ${APP_EXE} ${APP_PARAMS}

  # gather GPU power readings one last time
  nvidia-smi --query-gpu=index,timestamp,power.draw --format=csv &>> ${APP_RUN_PATH}/nvsmi-power-${SLURM_PROCID}.out

else

  ${APP_EXE} ${APP_PARAMS}

fi
