cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine ssqgh_dp
c		Vincent Lesur 13/06/2006
c
c       Modified 18.08.2011 to call ssqgh_d in place of ssqgh as
c       none diagonal covariance matrix element are not expected
c       for Large system (V.Lesur) 
c
c       Parallel interface for ssqgh.f
c
c       Called: ssqgh, MPI_ALLREDUCE
c
c       input:
c          shdeg          max SH degree value
c          nb             Number or base function to use
c          nlocpts        number of data+sampling points local to rank
c          nlocdatpts     number of data points assigned to rank
c          bc             Estimation of Base function coefficients
c          ddat           data vector
c          xyzf           result of forward modelling
c
c       output:
c          gj             gradient of the weighted sum of squares (nb)
c          dh             diagonal of the Hessian (nb)
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ssqgh_dp(shdeg, nb,
     >                      nlocpts, nlocdatpts,
     >                      bc, ddat, xyzf,
     >                      gj, dh)
c
        use cudafor, only : cudaDeviceSynchronize
        use kernels
c
        implicit none
c
        include 'mpif.h'
c
        integer shdeg, nb
        integer nlocpts, nlocdatpts
c
        real*8 bc(1:nb)        
        real*8 ddat(1:nlocpts)
        real*8 xyzf(1:nlocpts)
        real*8 gj(1:nb), dh(1:nb)
c
        integer ierr, istat
c
c
        call init_ssqgh_device_arrays(nb, nlocpts, bc, ddat, xyzf)
c
c
c  Offload data points
c          
c
#if defined(CUDA_KERNEL_LOOP)
c
        call ssqgh_dat_loop(shdeg, 1, nlocdatpts)
c      
#else       
c
        call ssqgh_dat_kernel <<< get_nblocks_dat(), get_nthreads() >>>
     >    (shdeg, 1, nlocdatpts, 0)
c
#endif
c    
        istat = cudaDeviceSynchronize()
c
#if defined(CUDA_DEBUG)
        call check_for_cuda_error('ssqgh_dat')
#endif
c
c
c  Offload sampling points
c  Workload is higher due to ssqgh_kernel_loop()/ssqgh_sam_kernel() calling XYZsph_bi0_sample()
c
c
#if defined(CUDA_KERNEL_LOOP)
c
        call ssqgh_sam_loop(shdeg, nlocdatpts+1, nlocpts)
c 
#else
c
        call ssqgh_sam_kernel <<< get_nblocks_sam(), get_nthreads() >>>
     >    (shdeg, nlocdatpts+1, nlocpts, nlocdatpts)
c
#endif
c    
        istat = cudaDeviceSynchronize()
c
#if defined(CUDA_DEBUG)
        call check_for_cuda_error('ssqgh_sam')
#endif
c
c
        call get_ssqgh_device_arrays(nb, gj, dh)
c        
c
        call MPI_ALLREDUCE(MPI_IN_PLACE, gj, nb,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
        call MPI_ALLREDUCE(MPI_IN_PLACE, dh, nb,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
c
        end subroutine ssqgh_dp