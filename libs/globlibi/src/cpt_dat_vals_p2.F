cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals_p2
c		Vincent Lesur 09/02/2005
c
c       Parallel interface for cpt_dat_vals.f
c
c       Called: cpt_dat_vals_p2
c
c       input:
c         nd            Space dimension
c         nlocpts       number of data+sampling points local to rank
c         nlocdatpts    number of data points local to rank
c         shdeg         max SH degree value
c         d2a           pre-computed array for mklf_F2()
c         ppos(nd+1,*)  point position in ndD
c         nb            Number of base functions
c         bc            base coefficients
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p2(nd, nlocpts, nlocdatpts, shdeg,
     >                             d2a, ppos, nb, bc, xyzf)
c
        implicit none
c
        integer :: nd, nb, nlocpts, nlocdatpts, shdeg
        real*8  :: d2a(0:shdeg), ppos(nd+1,nlocpts), bc(nb)
        real*8  :: wmam_fun, xyzf(nlocpts)
c
#ifdef OMP_OFFLOAD
        logical, save :: firstcall = .TRUE.
#endif
c
        integer :: i
c
        external wmam_fun
c 
c
c
#ifdef OMP_OFFLOAD
c
!$OMP TARGET DATA if(firstcall)
!$omp& map(to: nb, nd, nlocpts, nlocdatpts, shdeg)
!$omp& map(to: d2a(0:shdeg))
!$omp& map(to: ppos(1:nd+1,1:nlocpts))
        if (firstcall) then
          firstcall = .FALSE.
        endif
c
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!$omp& default(none)
!$omp& shared(nb, nd, nlocpts, nlocdatpts, shdeg)
!$omp& shared(d2a, ppos, bc, xyzf)
!$omp& map(to: bc(1:nb))
!$omp& map(from: xyzf(1:nlocpts))
!$omp& schedule(static)
        do i=1,nlocpts
          xyzf(i) = wmam_fun((i>nlocdatpts),
     >                       shdeg, nb, nd,
     >                       d2a, bc, ppos(1,i))
        enddo
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
!$OMP END TARGET DATA
c
#else
c
!$OMP PARALLEL DO
!$omp& default(shared)
!$omp& schedule(static)
        do i=1,nlocpts
          xyzf(i) = wmam_fun((i>nlocdatpts),
     >                       shdeg, nb, nd,
     >                       d2a, bc, ppos(1,i))
        enddo
!$OMP END PARALLEL DO
c
#endif
c
c
        return
        end
