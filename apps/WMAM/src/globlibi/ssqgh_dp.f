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
c          nd             space dimension
c          nlocpts        number of data+sampling points local to rank
c          nlocdatpts     number of data points assigned to rank
c          d2a            pre-computed array used by mk_lf_dlf()
c          bc             Estimation of Base function coefficients
c          ppos           data point position
c          cov            Covariance matrix in SLAP column format
c          jcov           integer arrays describing cov format
c          ddat           data vector
c          xyzf           result of forward modelling
c
c       output:
c          gj             gradient of the weighted sum of squares (nb)
c          dh             diagonal of the Hessian (nb)
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ssqgh_dp(shdeg, nb, nd,
     >                      nlocpts, nlocdatpts,
     >                      d2a, bc, ppos,
     >                      cov, jcov, ddat, xyzf,
     >                      gj, dh)
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
        real*8 cov(1:nlocpts)
        integer jcov(1:nlocpts+2)
        real*8 ddat(1:nlocpts)
        real*8 xyzf(1:nlocpts)
        real*8 gj(1:nb), dh(1:nb)
c
        real(8), allocatable, device :: d_d2a(:)
        real(8), allocatable, device :: d_bc(:)
        real(8), allocatable, device :: d_ppos(:,:)
        real(8), allocatable, device :: d_cov(:)
        integer, allocatable, device :: d_jcov(:)
        real(8), allocatable, device :: d_ddat(:)
        real(8), allocatable, device :: d_xyzf(:) 
        real(8), allocatable, device :: d_gj(:)
        real(8), allocatable, device :: d_dh(:)   
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
        allocate(d_cov(1:nlocpts))
        allocate(d_jcov(1:nlocpts+2))
        allocate(d_ddat(1:nlocpts))
        allocate(d_xyzf(1:nlocpts))
        allocate(d_gj(1:nb))
        allocate(d_dh(1:nb))
c
        d_d2a  = d2a
        d_bc   = bc
        d_ppos = ppos
        d_cov  = cov
        d_jcov = jcov
        d_xyzf = xyzf
        d_ddat = ddat
        d_gj(1:nb) = 0.0d0
        d_dh(1:nb) = 0.0d0
c
c
        nthreads = 128
        nblocks = nlocdatpts / nthreads
        if (MOD(nlocdatpts,nthreads) .gt. 0) then
          nblocks = nblocks + 1
        endif
c
        call ssqgh_dp_dat<<<nblocks,nthreads>>>
     >    (shdeg, nb, nd, nlocpts, nlocdatpts,
     >     d_d2a, d_ppos,
     >     d_cov, d_jcov, d_ddat, d_xyzf,
     >     d_gj, d_dh)
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
     >    (shdeg, nb, nd, nlocpts, nlocdatpts,
     >     d_d2a, d_bc, d_ppos,
     >     d_cov, d_jcov, d_ddat, d_xyzf,
     >     d_gj, d_dh)
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
        gj = d_gj
        dh = d_dh
c
c
        deallocate(d_d2a)
        deallocate(d_bc)
        deallocate(d_ppos)
        deallocate(d_cov)
        deallocate(d_jcov)
        deallocate(d_ddat)
        deallocate(d_xyzf)
        deallocate(d_gj)
        deallocate(d_dh)
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