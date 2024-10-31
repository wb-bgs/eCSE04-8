cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals_p
c		Vincent Lesur 09/02/2005
c
c       Parallel interface for cpt_dat_vals.f
c
c       Called: cpt_dat_vals_p
c
c       input:
c         nd            Space dimension
c         nlocpts       number of data+sampling points local to rank
c         nlocdatpts    number of data points local to rank
c         shdeg         max SH degree value
c         d2a           pre-computed array for mk_lf_dlf()
c         ppos(nd+1,*)  point position in nd
c         nb            Number of base functions
c         bc            base coefficients
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p(nd, nlocpts, nlocdatpts, shdeg,
     >                            d2a, ppos, nb, bc, xyzf)
c
        implicit none
c
        include 'mpif.h'
c
        integer :: nd, nb, nlocpts, nlocdatpts, shdeg
        real*8  :: d2a(0:shdeg), ppos(nd+1,nlocpts), bc(nb)
        real*8  :: xyzf(nlocpts)
c
        real*8, parameter :: RAG = 6371.2d0
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c        
        integer i
        real*8 p1, p2, ra
        real*8 bex, bey, bez
c
        real*8, allocatable :: dra(:)
        real*8, allocatable :: dalpha(:), dbeta(:)
        real*8, allocatable :: dlf(:), ddlf(:)
c
        real*8 XYZsph_bi0_fun
c
#if defined(OMP_OFFLOAD_CPTP)
        logical, save :: firstcall = .TRUE.
#endif
c
c   
        allocate(dra(shdeg))
        allocate(dalpha(0:shdeg), dbeta(0:shdeg))
        allocate(dlf(shdeg+1), ddlf(shdeg+1)) 
c
#if defined(OMP_OFFLOAD_CPTP)
!$OMP TARGET DATA if(firstcall)
!$omp& map(to: nb, nd)
!$omp& map(to: nlocpts, nlocdatpts, shdeg)
!$omp& map(to: d2a(0:shdeg), dra(1:shdeg))
!$omp& map(to: dalpha(0:shdeg), dbeta(0:shdeg))
!$omp& map(to: dlf(1:shdeg+1), ddlf(1:shdeg+1))
!$omp& map(to: ppos(1:nd+1,1:nlocpts))
        if (firstcall) then
          firstcall = .FALSE.
        endif
c
!$OMP TARGET DATA
!$omp& map(to: bc(1:nb))
#endif
c
c
#if defined(OMP_OFFLOAD_CPTP)
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!$omp& map(tofrom: xyzf(1:nlocdatpts))
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(dra, dalpha, dbeta)
!$omp& private(dlf, ddlf)
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
     >                              dalpha, dbeta,
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
!$omp& map(tofrom: xyzf(nlocdatpts+1:nlocpts))
#else
!$OMP PARALLEL DO
#endif
!$omp& default(shared)
!$omp& private(dra, dalpha, dbeta)
!$omp& private(dlf, ddlf)
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
     >                           dra, dalpha, dbeta,
     >                           dlf, ddlf, bc,
     >                           p1, p2, ra, 
     >                           bex, bey, bez)
c
          xyzf(i) = XYZsph_bi0_fun(shdeg, nb,
     >                             d2a, dra,
     >                             dalpha, dbeta,
     >                             dlf, ddlf, bc,
     >                             p1, p2, ra,
     >                             bex, bey, bez)
c
        enddo
#if defined(OMP_OFFLOAD_CPTP)
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
!$OMP END TARGET DATA
!$OMP END TARGET DATA
#else
!$OMP END PARALLEL DO
#endif
c
c
        deallocate(dra)
        deallocate(dalpha, dbeta)
        deallocate(dlf, ddlf)
c
        return
        end