ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          subroutine GC_STEP_P
c                           vincent Lesur 02.11.2010
c
c      18.08.2011 Triple-checked against DGC.F slatec subroutine
c        1- added the weight matrix W in calculation of zz(*) 
c        2- divide step by 2 because GJ=2*A'.W.R.W.[B-A.M] 
c
c      15.07.2011 (v.Lesur) update input list of cptstd
c      Modified:  10.3.2011 V.Lesur
c         1- sign of gj,ghj,ds changed 
c
c     Calculate the the length of a GC step, and make the associated 
c     Forward modelling
c
c     GC algorithm to solve A'.W.A.x=A'.W.b
c     x_new=x+step.ds
c     step = (gj.ghj)/(ds'.A'.W.A.ds)
c
c     Parallel version
c
c     input:
c         iunit         integer unit number for I/O
c         shdeg         max SH degree value
c         nb            Number or base function to use
c         nd            space dimension
c         npts          Total number of points (data + sampling) for all ranks
c         nlocpts       Total number of points for this rank
c         nlocdatpts    number of data points assigned to rank
c         d2a           pre-computed array used by mk_lf_dlf()
c         (d)dlf        pre-allocated arrays computed by mk_lf_dlf() and
c                       used within XYZsph_bi0
c         bc            Estimate of Base function coefficients
c         ppos          data point position in ndD
c         ddat          data values
c         cov           covariance matrix in SLAP Column format
c         jcov          Integer vector describing cov format
c         std           STD value for given BC
c         gj            gradient direction
c         ghj           preconditioned gradient direction
c         ds            descent direction
c
c       output:
c         stp           recommended step in direction ds(*)
c         std           STD value for given BC+stp*DS
c         xyzf          Forward modelling for given BC+stp*DS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine gc_step_p(iunit, shdeg, nb, nd, npts,
     >                       nlocpts, nlocdatpts,
     >                       d2a, dlf, ddlf,
     >                       bc, ppos, ddat,
     >                       cov, jcov,
     >                       std, gj, ghj,
     >                       ds, stp, xyzf)
c
        implicit none
c
        include 'mpif.h'
c
        integer iunit, shdeg, nb, nd
        integer npts, nlocpts, nlocdatpts
        real*8 d2a(0:shdeg)
        real*8 dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real*8 bc(1:nb)
        real*8 ppos(1:nd+1,1:nlocpts)
        real*8 ddat(1:nlocpts)
        real*8 cov(1:nlocpts)
        integer jcov(1:nlocpts+2)
        real*8 std, gj(1:nb), ghj(1:nb)
        real*8 ds(1:nb), stp, xyzf(1:nlocpts)
c
        integer i
        integer ierr, rank 
        real*8, allocatable :: bcn(:), zz(:)
        real*8 zzs
c        
c
c All: Defining parallel enviroment
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
c All: Calculate  sqrt(w).A.DS 
        allocate(zz(1:nlocpts))
        zz(1:nlocpts) = 0.0d0
        call cpt_dat_vals_p(shdeg, nb, nd, nlocpts, nlocdatpts,
     >                      d2a, dlf, ddlf, ds, ppos, zz)
c
        do i = 1,nlocpts
          zz(i) = zz(i)/dsqrt(cov(jcov(i)))
          zz(i) = zz(i)**2
        enddo
        zzs = SUM(zz)
c
c  All: calculate step
        stp = dot_product(gj(1:nb), ghj(1:nb))

        call MPI_ALLREDUCE(MPI_IN_PLACE, zzs, 1,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)

        stp = stp/zzs
        stp = stp/2.d0
        deallocate(zz)
c
        if (rank.eq.0) then
            write(iunit,*) 'GC_STEP :',stp
        endif
c
c ALL: Estimate the new set of parameter for a step stp in direction ds
        allocate(bcn(1:nb))
        bcn(1:nb) = bc(1:nb) + stp*ds(1:nb)
c
c ALL: Do the forward modelling
        xyzf(1:nlocpts) = 0.0d0
        call cpt_dat_vals_p(shdeg, nb, nd, nlocpts, nlocdatpts,
     >                      d2a, dlf, ddlf, bcn, ppos, xyzf)
c
        call cptstd_dp(npts, nlocpts,
     >                 cov, jcov, ddat,
     >                 xyzf, std)
c
        deallocate(bcn)
c
        return
        end