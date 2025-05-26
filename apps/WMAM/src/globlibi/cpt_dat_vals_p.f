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
c         (d)dlf        pre-allocated arrays computed by mk_lf_dlf() and
c                       used within XYZsph_bi0
c         bc            base coefficients
c         ppos(nd+1,*)  point position in nd
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p(shdeg, nb, nd, nlocpts, nlocdatpts,
     >                            d2a, dra, dlf, ddlf, bc, ppos, xyzf)
c
        implicit none
c
        integer shdeg, nb, nd, nlocpts, nlocdatpts
        real*8  d2a(0:shdeg), dra(1:shdeg)
        real*8  dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real*8  bc(1:nb), ppos(nd+1,nlocpts)
        real*8  xyzf(1:nlocpts)
c
        real*8, parameter :: RAG = 6371.2d0
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c        
        integer i
        real*8 p1, p2, ra
        real*8 bex, bey, bez
c
        real*8 XYZsph_bi0_fun
c
c
#if defined(OMP_OFFLOAD_CPTP)
!$OMP TARGET DATA
!$omp& map(to: bc(1:nb))
!$omp& map(tofrom: xyzf(1:nlocpts))
#endif
c
c
#if defined(OMP_OFFLOAD_CPTP)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(dra, dlf, ddlf)
!$omp& private(p1, p2, ra)
!$omp& private(bex, bey, bez)
#if defined(OMP_OFFLOAD_CPTP)
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
c
           xyzf(i) = XYZsph_bi0_fun(shdeg, nb, d2a,
     >                              dra, dlf, ddlf,
     >                              bc, p1, p2, ra,
     >                              bex, bey, bez)
c
        enddo
#if defined(OMP_OFFLOAD_CPTP)
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP END PARALLEL DO
#endif
c
c
#if defined(OMP_OFFLOAD_CPTP)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(dra, dlf, ddlf)
!$omp& private(p1, p2, ra)
!$omp& private(bex, bey, bez)
#if defined(OMP_OFFLOAD_CPTP)
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
          call XYZsph_bi0_sample(shdeg, nb, d2a,
     >                           dra, dlf, ddlf,
     >                           bc, p1, p2, ra, 
     >                           bex, bey, bez)
c
          xyzf(i) = XYZsph_bi0_fun(shdeg, nb, d2a,
     >                             dra, dlf, ddlf,
     >                             bc, p1, p2, ra,
     >                             bex, bey, bez)
c
        enddo
#if defined(OMP_OFFLOAD_CPTP)
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
c
!$OMP END TARGET DATA
#else
!$OMP END PARALLEL DO
#endif
c
c
        return
        end