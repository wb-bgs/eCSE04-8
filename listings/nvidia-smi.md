# Reading NVIDIA GPU power

The `nvidia-smi` program is called from each GPU node assigned to the job. Something like the following script is called from the Slurm `srun` command.


```bash

...

if ! ((${SLURM_LOCALID} % ${MPI_RANKS_PER_NODE})); then

  nvidia-smi --query-gpu=index,timestamp,power.draw.instant --loop-ms=1000 \
      --format=csv &> nvsmi-power-${SLURM_PROCID}.out &
  
  ${APP_EXE} ${APP_PARAMS}

  nvidia-smi --query-gpu=index,timestamp,power.draw.instant \
      --format=csv &>> nvsmi-power-${SLURM_PROCID}.out

else

  ${APP_EXE} ${APP_PARAMS}

endif

...

```
