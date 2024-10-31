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
c          nd             space dimension
c          nlocpts        number of data+sampling points local to rank
c          nlocdatpts     number of data points assigned to rank
c          shdeg          max SH degree value
c          d2a            pre-computed array for mk_lf_dlf()
c          ppos           data point position in ndD
c          nb             Number or base function to use
c          bc             Estimation of Base function coefficients
c          jcov           integer arrays describing cov format
c          cov            Covariance matrix in SLAP column format
c          ddat           data vector
c          xyzf           result of forward modelling
c
c       output:
c          gj             gradient of the weighted sum of squares (nb)
c          dh             diagonal of the Hessian (nb)
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ssqgh_dp(nd, nlocpts, nlocdatpts, shdeg,
     >                      d2a, ppos, nb, bc,
     >                      jcov, cov, ddat, xyzf,
     >                      gj_map_len, gj_map,
     >                      gj, dh)
c
        implicit none
c
        include 'mpif.h'
c
        integer nd, nlocpts, nlocdatpts
        integer shdeg, nb, jcov(nlocpts+2)
        real*8 d2a(0:shdeg), ddat(nlocpts)
        real*8 xyzf(nlocpts), cov(nlocpts)
        real*8 ppos(nd+1,nlocpts), bc(nb)
        integer gj_map_len 
        integer gj_map(gj_map_len)
        real*8 gj(nb), dh(nb)
c
        real*8, parameter :: RAG = 6371.2d0
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
        integer i, j, nu, ierr
        real*8 p1, p2, ra
        real*8 bex, bey, bez
c
        real*8, allocatable :: dra(:)
        real*8, allocatable :: dalpha(:), dbeta(:)
        real*8, allocatable :: dlf(:), ddlf(:)
c
        real*8 dw_dh, dw_gj
c
        real*8, allocatable :: gj2(:), dh2(:)
c 
#if defined(OMP_OFFLOAD_SSQGH)
        logical, save :: firstcall = .TRUE.
#endif
c
c
        allocate(dra(shdeg))
        allocate(dalpha(0:shdeg), dbeta(0:shdeg))
        allocate(dlf(shdeg+1), ddlf(shdeg+1))
c 
        allocate(gj2(nb))
        allocate(dh2(nb))
c
        gj2(1:nb) = 0.0d0
        dh2(1:nb) = 0.0d0
c
c       
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP TARGET DATA if (firstcall)
#if !defined(OMP_OFFLOAD_CPTP)
!$omp& map(to: nb, nd)
!$omp& map(to: nlocpts, nlocdatpts, shdeg)
!$omp& map(to: d2a(0:shdeg), dra(1:shdeg))
!$omp& map(to: dalpha(0:shdeg), dbeta(0:shdeg))
!$omp& map(to: dlf(1:shdeg+1), ddlf(1:shdeg+1))
!$omp& map(to: ppos(1:nd+1,1:nlocpts))
#endif
!$omp& map(to: cov(1:nlocpts), jcov(1:nlocpts+2))
        if (firstcall) then
          firstcall = .FALSE.
        endif
c
!$OMP TARGET DATA
!$omp& map(to: ddat(1:nlocpts))
!$omp& map(to: xyzf(1:nlocpts))
!$omp& map(to: bc(1:nb))
!$omp& map(tofrom: gj2(1:nb), dh2(1:nb))
#endif
c
c
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(dra, dalpha, dbeta)
!$omp& private(dlf, ddlf)
!$omp& private(p1, p2, ra)
!$omp& private(bex, bey, bez)
!$omp& private(dw_dh, dw_gj)
#if defined(OMP_OFFLOAD_SSQGH)
!$omp& dist_schedule(static)
#else
!$omp& schedule(static)
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
     >                        dra, dalpha, dbeta,
     >                        dlf, ddlf,
     >                        p1, p2, ra,
     >                        bex, bey, bez,
     >                        dw_gj, dw_dh,
     >                        gj2, dh2)
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
!$omp& private(dra, dalpha, dbeta)
!$omp& private(dlf, ddlf)
!$omp& private(p1, p2, ra)
!$omp& private(bex, bey, bez)
!$omp& private(dw_dh, dw_gj)
#if defined(OMP_OFFLOAD_SSQGH)
!$omp& dist_schedule(static)
#else
!$omp& schedule(static)
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
          call XYZsph_bi0_sample(shdeg, nb,
     >                           d2a, dra,
     >                           dalpha, dbeta,
     >                           dlf, ddlf, bc,
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
     >                        dra, dalpha, dbeta,
     >                        dlf, ddlf,
     >                        p1, p2, ra,
     >                        bex, bey, bez,
     >                        dw_gj, dw_dh,
     >                        gj2, dh2)
c
        enddo
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
!$OMP END TARGET DATA
!$OMP END TARGET DATA
#else
!$OMP END PARALLEL DO
#endif
c
c
        call MPI_ALLREDUCE(MPI_IN_PLACE, gj2, nb,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
        call MPI_ALLREDUCE(MPI_IN_PLACE, dh2, nb,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
c
c  Rearrange gj/dh coefficient terms
        do i = 1,shdeg
          nu = gj_map(i)
          gj(nu) = gj2(i)
          dh(nu) = dh2(i)
        enddo
c
        j = shdeg+1
        do i = shdeg+1,gj_map_len
          nu = gj_map(i)
          gj(nu) = gj2(j)
          gj(nu+1) = gj2(j+1)
          dh(nu) = dh2(j)
          dh(nu+1) = dh2(j+1)
          nu = nu+2
          j = j+2
        enddo
c
        deallocate(gj2, dh2)
c
        deallocate(dra)
        deallocate(dalpha, dbeta)
        deallocate(dlf, ddlf)
c
        return
        end