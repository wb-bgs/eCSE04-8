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
c         nd            Space dimension
c         nlocpts       number of data+sampling points local to rank
c         nlocdatpts    number of data points local to rank
c         d2a           pre-computed array used by mk_lf_dlf()
c         bc            base coefficients
c         ppos(nd+1,*)  point position in nd
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p(shdeg, nb, nd, nlocpts, nlocdatpts,
     >                            d2a, bc, ppos, xyzf)
c
        use cudafor
        use kernels
c
        implicit none
c
        include 'mpif.h'
c
        integer shdeg, nb, nd
        integer nlocpts, nlocdatpts
c
        real*8 d2a(0:shdeg)
        real*8 bc(1:nb)
        real*8 ppos(1:nd+1,1:nlocpts)
        real*8 xyzf(1:nlocpts)
c
        real(8), allocatable, device :: d_d2a(:)
        real(8), allocatable, device :: d_bc(:)
        real(8), allocatable, device :: d_ppos(:,:)
        real(8), allocatable, device :: d_xyzf(:)
c
        integer nthreads, nblocks, nlocsampts
        integer rank, istat, ierr     
c
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
c
        allocate(d_d2a(0:shdeg))
        allocate(d_bc(1:nb))
        allocate(d_ppos(1:nd+1,1:nlocpts))
        allocate(d_xyzf(1:nlocpts))
c
        d_d2a  = d2a
        d_bc   = bc
        d_ppos = ppos
c
c
        nthreads = 128
        nblocks = nlocdatpts / nthreads
        if (MOD(nlocdatpts,nthreads) .gt. 0) then
          nblocks = nblocks + 1
        endif
c
        call cpt_dat_vals_p_dat<<<nblocks,nthreads>>>
     >    (shdeg, nb, nd, nlocpts, nlocdatpts,
     >     d_d2a, d_bc, d_ppos, d_xyzf)
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
     >    (shdeg, nb, nd, nlocpts, nlocdatpts,
     >     d_d2a, d_bc, d_ppos, d_xyzf)

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
        xyzf = d_xyzf
c
c
        deallocate(d_d2a)
        deallocate(d_bc)
        deallocate(d_ppos)
        deallocate(d_xyzf)
c
c
        end subroutine cpt_dat_vals_p