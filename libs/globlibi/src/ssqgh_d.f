cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine ssqgh_d
c		Vincent Lesur  08/06/2006
c
c      Modified 14.07.2011 & 17.08.2011 (V. Lesur)
c       1- input modified in concoct subs
c       3- for diagonal covariance matrix only
c      Modified 10.3.2011 (V. Lesur)
c       1- sign of ddif (and therefore GJ) changed
c       2- Scaling by npt of Gj and Hj removed
c
c       computes the Gradient of the weighted sum of the squared 
c       differences between data and model
c       Is also calulated the associated diagonal of the Hessian matrix
c
c       Should be used with F90 ot later versions
c
c       Compute the normal equations of the problem calling recursive 
c       subroutines
c
c       input:
c          npmax          number max of data points handled together
c          nlocpts        number of data+sampling points local to rank
c          nlocdatpts     number of data points assigned to rank
c          shdeg          max SH degree value
c          d2a            pre-computed array for mklf_F2()
c          nd             space dimension
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
c          hj             diagonal of the Hessian (nb)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ssqgh_d(npmax, nlocpts, nlocdatpts,
     >                     shdeg, d2a, nd, ppos, nb,
     >                     bc, jcov, cov,
     >                     ddat, xyzf, gj, hj)
c
        implicit none
c
        integer npmax, nlocpts, nd, nb, nlocdatpts, shdeg
        integer jcov(nlocpts+2)
        real*8 d2a(0:shdeg), ppos(nd+1,nlocpts), bc(nb)
        real*8 cov(nlocpts), ddat(nlocpts)
        real*8 xyzf(nlocpts), gj(nb), hj(nb)
c
#ifdef OMP_OFFLOAD
        logical, save :: firstcall = .TRUE.
#endif
c
        integer nchk, nrem
        integer ip, ip2, i, j
        real*8, allocatable :: dwgh(:),ddif(:)
        real*8, allocatable :: aa(:,:)
c
c
        gj(1:nb)=0.0d0
        hj(1:nb)=0.0d0
c
        nchk = nlocpts / npmax
        nrem = MOD(nlocpts, npmax)
c
        allocate(dwgh(npmax),ddif(npmax))
        allocate(aa(nb,npmax))
c
c
#ifdef OMP_OFFLOAD
c
c
!$OMP TARGET DATA if (firstcall)
!$omp& map(to: nb, nd, npmax)
!$omp& map(to: nlocpts, nlocdatpts, shdeg)
!$omp& map(to: d2a(0:shdeg))
!$omp& map(to: ppos(1:nd+1,1:nlocpts))
!$omp& map(to: cov(1:nlocpts), jcov(1:nlocpts+2))
!$omp& map(alloc: dwgh(1:npmax), ddif(1:npmax))
!$omp& map(alloc: aa(1:nb,1:npmax))
        if (firstcall) then
          firstcall = .FALSE.
        endif
c
!$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO
!$omp& default(none)
!$omp& shared(nb, nd, npmax, nchk)
!$omp& shared(nlocpts, nlocdatpts, shdeg)
!$omp& shared(d2a, ppos, cov, jcov, bc)
!$omp& shared(ddat, xyzf, aa)
!$omp& private(ip, ip2, j)
!$omp& private(dwgh, ddif)
!$omp& map(to: bc(1:nb))
!$omp& map(to: ddat(1:nlocpts), xyzf(1:nlocpts))
!$omp& map(tofrom: gj(1:nb), hj(1:nb))
!$omp& schedule(static)
!$omp& reduction(+:gj,hj)
        do i = 1,nchk
          ip = 1 + (i-1)*npmax

          ip2 = ip
          do j = 1,npmax
            ddif(j) = ddat(ip2)-xyzf(ip2)
            dwgh(j) = 1.d0/cov(jcov(ip2))
            ip2 = ip2+1
          enddo
c        
c  calculate the equations of condition          
          call mkArows((ip>nlocdatpts),
     >                  shdeg, nb, nd, npmax, npmax,
     >                  d2a, bc, ppos(1,ip), aa)
c
c  update the G matrix and B vector
          call concoct_GJ(nb,npmax,npmax,dwgh,aa,ddif,gj)
          call concoct_HJ(nb,npmax,npmax,dwgh,aa,hj)
        enddo
!$OMP END TARGET TEAMS DISTRIBUTE PARALLEL DO
!$OMP END TARGET DATA
c
c
#else
c
c
!$OMP PARALLEL DO
!$omp& default(shared)
!$omp& private(ip,ip2,j)
!$omp& private(dwgh,ddif,aa)
!$omp& schedule(static)
!$omp& reduction(+:gj,hj)
        do i = 1,nchk
          ip = 1 + (i-1)*npmax

          ip2 = ip
          do j = 1,npmax
            ddif(j) = ddat(ip2)-xyzf(ip2)
            dwgh(j) = 1.d0/cov(jcov(ip2))
            ip2 = ip2+1
          enddo
c        
c  calculate the equations of condition          
          call mkArows((ip>nlocdatpts),
     >                  shdeg, nb, nd, npmax, npmax,
     >                  d2a, bc, ppos(1,ip), aa)
c
c  update the G matrix and B vector
          call concoct_GJ(nb,npmax,npmax,dwgh,aa,ddif,gj)
          call concoct_HJ(nb,npmax,npmax,dwgh,aa,hj)
        enddo
!$OMP END PARALLEL DO
c
c
#endif
c
c
c  do the last iteration
        ip = 1 + nchk*npmax
        ip2 = ip
        do j = 1,nrem
          ddif(j) = ddat(ip2)-xyzf(ip2)
          dwgh(j) = 1.d0/cov(jcov(ip2))
          ip2 = ip2+1
        enddo
c
c  calculate the equations of condition
        call mkArows((ip>nlocdatpts),
     >               shdeg, nb, nd, nrem, npmax,
     >               d2a, bc, ppos(1,ip), aa)
c
c  update the G matrix and B vector
        call concoct_GJ(nb,nrem,npmax,dwgh,aa,ddif,gj)
        call concoct_HJ(nb,nrem,npmax,dwgh,aa,hj)
c
c
        deallocate(dwgh,ddif,aa)
c
c
        return
        end
