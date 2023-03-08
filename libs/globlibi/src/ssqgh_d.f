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
c          ipg            where to start in data file!
c          nd             space dimension
c          ppos           data point position in ndD
c          nb             Number or base function to use
c          fun_mf         misfit function (like l2_norm.f)
c          sub_base       the "Base functions" subroutine to use
c          bc             Estimation of Base function coefficients
c          jcov           integer arrays describing cov format
c          cov            Covariance matrix in SLAP column format
c          ddat           data vector
c          ntv(2)         basis number for each point
c          ntn(2)         point count for each basis number
c          xyzf           result of forward modelling
c
c       output:
c          gj             gradient of the weighted sum of squares (nb)
c          hj             diagonal of the Hessian (nb)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ssqgh_d(npmax, nlocpts, ipg,
     >                     nd, ppos, nb,
     >                     fun_mf, sub_base, bc, jcov, cov,
     >                     ddat, xyzf, gj, hj)
c
        implicit none
c
        integer ip,nlocpts
        integer np,npmax,nd,nb,jcov(*)
        integer i,ipg,ipl
        integer, allocatable :: ntval(:)
        real*8 ddat(*),xyzf(*),cov(*),ppos(nd+1,*),bc(*)
        real*8 gj(*),hj(*)
c       real*8, allocatable :: vmf(:)
        real*8, allocatable :: dwgh(:),ddif(:)
        real*8, allocatable :: aa(:,:)
c
        real*8 fun_mf
        external fun_mf,sub_base
c
c  ipg : ip global
c  ipl : ip local
c
c       allocate(vmf(npmax))
        allocate(dwgh(npmax),ddif(npmax))
        allocate(aa(nb,npmax))
        allocate(ntval(npmax))
c
        gj(1:nb)=0.0d0
        hj(1:nb)=0.0d0
        ntval(1:npmax)=0
c
c
        ip=1
        do while (ip.le.nlocpts) 
          ipl=ipg+ip-1
          np=min0(nlocpts-ip+1,npmax)

          do i=1,np
            ntval(i)=ip-1+i
          enddo
c        
c
c  calculate the equations of condition
          call mkArows(np,ntval,nd,nb,ppos(1,ipl),sub_base,bc,aa)
c
c  calculate the delta data
          do i=1,np
            ddif(i)=ddat(ipl+i-1)-xyzf(ipl+i-1)
          enddo
c
c  calculate the inverse covariance matrix
          do i=1,np
            dwgh(i)=1.d0/cov(jcov(ipl+i-1))
          enddo
c
c  calculate msft vector
c         do i=1,np
c            vmf(i)=ddif(i)*dsqrt(dwgh(i))
c         enddo

c
c  update the G matrix and B vector
          call concoct_GJ(fun_mf,nb,np,dwgh,aa,ddif,gj)
          call concoct_HJ(fun_mf,nb,np,dwgh,aa,hj)
c
          ip = ip + np
        enddo
c
        deallocate(dwgh,ddif,aa)
        deallocate(ntval)
c       deallocate(vmf)
c
        return
        end