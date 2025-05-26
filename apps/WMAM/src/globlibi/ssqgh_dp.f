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
c          (d)dlf         pre-allocated arrays computed by mk_lf_dlf() and
c                         used within XYZsph_bi0
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
     >                      d2a, dra, dlf, ddlf, bc, ppos,
     >                      cov, jcov, ddat, xyzf,
     >                      gj, dh)
c
        implicit none
c
        include 'mpif.h'
c
        integer shdeg, nb, nd, nlocpts, nlocdatpts
        real*8 d2a(0:shdeg), dra(1:shdeg)
        real*8 dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real*8 bc(1:nb)
        real*8 ppos(1:nd+1,1:nlocpts)
        real*8 cov(1:nlocpts)
        integer jcov(1:nlocpts+2)
        real*8 ddat(1:nlocpts)
        real*8 xyzf(1:nlocpts)
        real*8 gj(1:nb), dh(1:nb)
c
        real*8, parameter :: RAG = 6371.2d0
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
        integer i, ierr
        real*8 p1, p2, ra
        real*8 bex, bey, bez
c
        real*8 dw_dh, dw_gj
c
c
        gj(1:nb) = 0.0d0
        dh(1:nb) = 0.0d0
c
c       
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP TARGET DATA
!$omp& map(to: ddat(1:nlocpts))
!$omp& map(to: xyzf(1:nlocpts))
!$omp& map(to: bc(1:nb))
!$omp& map(tofrom: gj(1:nb), dh(1:nb))
#endif
c
c
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(dra, dlf, ddlf)
!$omp& private(p1, p2, ra)
!$omp& private(bex, bey, bez)
!$omp& private(dw_dh, dw_gj)
#if defined(OMP_OFFLOAD_SSQGH)
!$omp& dist_schedule(static)
#else
!$omp& schedule(static)
!$omp& reduction(+:gj,dh)
#endif
        do i = 1,nlocdatpts
c       
          p1 = ppos(1,i)*D2R
          p2 = ppos(2,i)*D2R
          ra = RAG / ppos(3,i)
c
          bex = ppos(5,i)
          bey = ppos(6,i)
          bez = ppos(7,i)

c  calculate the equations of condition   
c  and update the G matrix and B vector 
c
          dw_dh = 2.d0*(1.d0/cov(jcov(i)))
          dw_gj = dw_dh*(ddat(i)-xyzf(i))
c      
          call XYZsph_bi0_sub(shdeg, nb, d2a,
     >                        dra, dlf, ddlf,
     >                        p1, p2, ra,
     >                        bex, bey, bez,
     >                        dw_gj, dw_dh,
     >                        gj, dh)
c
        enddo
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP END PARALLEL DO
#endif
c
c
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(dra, dlf, ddlf)
!$omp& private(p1, p2, ra)
!$omp& private(bex, bey, bez)
!$omp& private(dw_dh, dw_gj)
#if defined(OMP_OFFLOAD_SSQGH)
!$omp& dist_schedule(static)
#else
!$omp& schedule(static)
!$omp& reduction(+:gj,dh)
#endif
        do i = nlocdatpts+1,nlocpts
c
          p1 = ppos(1,i)*D2R
          p2 = ppos(2,i)*D2R
          ra = RAG / ppos(3,i)
c
          bex = ppos(5,i)
          bey = ppos(6,i)
          bez = ppos(7,i)
c
          call XYZsph_bi0_sample(shdeg, nb, d2a,
     >                           dra, dlf, ddlf, bc,
     >                           p1, p2, ra, 
     >                           bex, bey, bez)
c        
c  calculate the equations of condition   
c  and update the G matrix and B vector 
c
          dw_dh = 2.d0*(1.d0/cov(jcov(i)))
          dw_gj = dw_dh*(ddat(i)-xyzf(i))
c      
          call XYZsph_bi0_sub(shdeg, nb, d2a,
     >                        dra, dlf, ddlf,
     >                        p1, p2, ra,
     >                        bex, bey, bez,
     >                        dw_gj, dw_dh,
     >                        gj, dh)
c
        enddo
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
c
!$OMP END TARGET DATA
#else
!$OMP END PARALLEL DO
#endif
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
        return
        end
