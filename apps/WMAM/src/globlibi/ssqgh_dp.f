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
        use cudafor
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
        integer nthreads, nblocks, nlocsampts
        integer rank, istat, ierr 
c
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
c
        call init_ssqgh_device_arrays(nb, nlocpts,
     >                                bc, ddat, xyzf)
c
c
        nthreads = 128
        nblocks = nlocdatpts / nthreads
        if (MOD(nlocdatpts,nthreads) .gt. 0) then
          nblocks = nblocks + 1
        endif
c
        call ssqgh_dp_dat<<<nblocks,nthreads>>>
     >    (shdeg, nlocpts, nlocdatpts)
c
        istat = cudaDeviceSynchronize()
c
#if defined(CUDA_DEBUG)
        ierr = cudaGetLastError()
        if (ierr .gt. 0) then
          write(*,*) rank,
     >        ': Error, ssqgh_dp_dat kernel failure: ',
     >        ierr, ', ', cudaGetErrorString(ierr)
          stop
        endif
#endif
c
c
        nthreads = 128
        nlocsampts = nlocpts - nlocdatpts
        nblocks = nlocsampts / nthreads
        if (MOD(nlocsampts,nthreads) .gt. 0) then
          nblocks = nblocks + 1
        endif
c
        call ssqgh_dp_smp<<<nblocks,nthreads>>>
     >    (shdeg, nlocpts, nlocdatpts)
c
        istat = cudaDeviceSynchronize()
c
#if defined(CUDA_DEBUG)
        ierr = cudaGetLastError()
        if (ierr .gt. 0) then
          write(*,*) rank,
     >        'Error, ssqgh_dp_smp kernel failure: ',
     >        ierr, ', ', cudaGetErrorString(ierr)
          stop
        endif
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