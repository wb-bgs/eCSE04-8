cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals_p
c		Vincent Lesur 09/02/2005
c
c       Parallel interface for cpt_dat_vals.f
c
c       Called: cpt_dat_vals_p
c
c       input:
c         shdeg         max SH degree value
c         nb            Number of base functions
c         nlocpts       number of data+sampling points local to rank
c         nlocdatpts    number of data points local to rank
c         bc            base coefficients
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p(shdeg, nb, nlocpts, nlocdatpts,
     >                            bc, xyzf)
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
        real*8 xyzf(1:nlocpts)
c
        integer nthreads, nblocks, nlocsampts
        integer rank, istat, ierr     
c
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
c
        call init_cpt_device_arrays(nb, bc)
c
c
        nthreads = 128
        nblocks = nlocdatpts / nthreads
        if (MOD(nlocdatpts,nthreads) .gt. 0) then
          nblocks = nblocks + 1
        endif
c
        call cpt_dat_vals_p_dat<<<nblocks,nthreads>>>
     >    (shdeg, nlocpts, nlocdatpts)
c
        istat = cudaDeviceSynchronize()
c
#if defined(CUDA_DEBUG)
        ierr = cudaGetLastError()
        if (ierr .gt. 0) then
          write(*,*) rank,
     >        ': Error, cpt_dat_vals_p_dat kernel failure: ',
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
        call cpt_dat_vals_p_smp<<<nblocks,nthreads>>>
     >    (shdeg, nlocpts, nlocdatpts)

        istat = cudaDeviceSynchronize()
c
#if defined(CUDA_DEBUG)
        ierr = cudaGetLastError()
        if (ierr .gt. 0) then
          write(*,*) rank,
     >        ': Error, cpt_dat_vals_p_smp kernel failure: ',
     >        ierr, ', ', cudaGetErrorString(ierr)
          stop
        endif
#endif
c
c
        call get_cpt_device_arrays(nlocpts, xyzf)
c
c
        end subroutine cpt_dat_vals_p