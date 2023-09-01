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
c          nd             space dimension
c          ppos           data point position in ndD
c          nb             Number or base function to use
c          fun_mf         misfit function (like l2_norm.f)
c          sub_base       the "Base functions" subroutine to use
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
        subroutine ssqgh_d(npmax, nlocpts,
     >                     nd, ppos, nb,
     >                     fun_mf, sub_base, bc, jcov, cov,
     >                     ddat, xyzf, gj, hj)
c
        implicit none
c
        integer rank, npmax, nlocpts, nd, nb
        integer jcov(*)
        real*8 ppos(nd+1,*), bc(*), cov(*), ddat(*)
        real*8 xyzf(*), gj(nb), hj(nb)
c        
        real*8 fun_mf
        external fun_mf, sub_base
c
        integer nchk, nrem
        integer np, ip, ip2, i, j
        real*8, allocatable :: dwgh(:),ddif(:)
        real*8, allocatable :: aa(:,:)
c
c
        gj(1:nb)=0.0d0
        hj(1:nb)=0.0d0
c
        nchk = nlocpts / npmax
        nrem = MOD(nlocpts, npmax)
        np = npmax
c
c
!$OMP PARALLEL
!$OMP& DEFAULT(NONE)
!$OMP& SHARED(nchk,nrem,npmax,nd,nb)
!$OMP& SHARED(ddat,xyzf,cov,jcov,ppos,bc)
!$OMP& SHARED(gj,hj)
!$OMP& FIRSTPRIVATE(np)
!$OMP& PRIVATE(i,ip,ip2,j,dwgh,ddif,aa)
        allocate(dwgh(npmax),ddif(npmax))
        allocate(aa(nb,npmax))
c
!$OMP DO
!$OMP& SCHEDULE(STATIC)
!$OMP& REDUCTION(+:gj(1:nb),hj(1:nb))
        do i = 1,nchk+1
          if (i .gt. nchk) then
            np = nrem
          endif
          ip = 1 + (i-1)*npmax

          ip2 = ip
          do j = 1,np
            ddif(j) = ddat(ip2)-xyzf(ip2)
            dwgh(j) = 1.d0/cov(jcov(ip2))
            ip2 = ip2+1
          enddo
c        
c  calculate the equations of condition          
          call mkArows(np,ip,nd,nb,ppos(1,ip),sub_base,bc,aa)
c
c  update the G matrix and B vector
          call concoct_GJ(fun_mf,nb,np,dwgh,aa,ddif,gj)
          call concoct_HJ(fun_mf,nb,np,dwgh,aa,hj)
        enddo
!$OMP END DO
c
        deallocate(dwgh,ddif,aa)
!$OMP END PARALLEL
c
c
        return
        end