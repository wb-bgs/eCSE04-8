# Loop Offload

The array iterates over the `nlocpts` data points assigned to the MPI rank.
These points are made up of input puts (`nlocinppts`) and sampled points (`nlocsampts`).
The MPI ranks have equal numbers of input and sampled points.

```fortran
subroutine cpt_dat_vals_p(...)

  ...
  
  do i = 1,nlocpts

    ...

  enddo

  ...

end subroutine cpt_dat_vals_p
```


For offloading to GPU, the loop is split such that there are distinct kernels for
the input and sampled points. This ensures that the threads of each kernel are performing
the same amount of work.

```fortran
subroutine cpt_dat_vals_p(...)

  ...
  
  ! Iterate over the input points

  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
  ...
  do i = 1,nlocinppts

    ...

  enddo
    
  ...

  ! Iterate over the sampled points

  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
  ...
  do i = nlocinppts+1, nlocpts

    ...

  enddo
  
  ...

end subroutine cpt_dat_vals_p
```