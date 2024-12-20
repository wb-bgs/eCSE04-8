#!/bin/bash

MPI_RANKS_PER_NODE=$1
APP_EXE=$2
APP_PARAMS=$3

NCU_PARAMS="--launch-skip 2 --launch-count 10 --kernel-name regex:\"cpt_dat_vals_p|ssqgh_dp\" --section MemoryWorkloadAnalysis --section SourceCounters --section SpeedOfLight --section SpeedOfLight_RooflineChart --export wmam-baseline-${SLURM_PROCID}"

if ! ((${SLURM_LOCALID} % ${MPI_RANKS_PER_NODE})); then
  ncu ${NCU_PARAMS} ${APP_EXE} ${APP_PARAMS}
else
  ${APP_EXE} ${APP_PARAMS}
fi
