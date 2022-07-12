cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cptAtWAds
c		Vincent Lesur  18/08/2011
c
c      derived from ssqgh_d.f, for diagonal COV matrix  only
c       Should be used with F90 ot later versions
c
c       input:
c          npmax          number max of data points handled together
c          npt            total number of data points
c          ipg            where to start in data file!
c          nd             space dimension
c          ppos           data point position in ndD
c          nb             Number or base function to use
c          fun_mf         misfit function (like l2_norm.f)
c          sub_base       the "Base functions" subroutine to use
c          bc             Estimation of Base function coefficients
c          ds             current descent direction
c          icov/jcov      integer arrays describing cov format
c          cov            Covariance matrix in SLAP column format
c          ddat           data vector
c          nt             vector indicating data type 
c          xyzf           result of forward modelling
c
c       output:
c          zz             Vector A^t.W.A.ds (nb)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cptAtWAds(npmax, npt, ipg, nd, ppos, nb,
     >                       fun_mf, sub_base, bc, ds,
     >                       icov, jcov, cov, ddat, nt,
     >                       xyzf, zz)
c
        implicit none
c
        integer ip,npt,np,npmax,nd,nb,icov(*),jcov(*)
        integer i,ipg,ipl,nt(*)
        real*8 ddat(*),xyzf(*),cov(*),ppos(nd+1,*),bc(*),ds(*)
        real*8 zz(*)
        real*8, allocatable :: dwgh(:),vmf(:)
        real*8, allocatable :: ddif(:),aa(:)
c
        real*8 fun_mf
        external fun_mf,sub_base
c
c  ipg : ip global
c  ipl : ip local
c
        allocate(dwgh(1:npmax),vmf(1:npmax))
        allocate(ddif(1:npmax),aa(1:npmax*nb))
c
        zz(1:nb)=0.0d0
c
        ip=1
        do while (ip.le.npt) 
          ipl=ipg+ip-1
          np=min0(npt-ip+1,npmax)
c
c Calculate the equations of condition
          call mkArows(np,nt(ipl),nd,nb,ppos(1,ipl),sub_base,bc,aa)
c
c  calculate the delta data
          do i=1,np
            ddif(i)=ddat(ipl+i-1)-xyzf(ipl+i-1)
          enddo
c
c  Calculate the inverse covariance matrix
          do i=1,np
            dwgh(i)=1.d0/cov(jcov(ipl+i-1))
          enddo
c
c Calculate msft vector
          do i=1,np
            vmf(i)=ddif(i)*dsqrt(dwgh(i))
          enddo
c
c Update the G matrix and B vector
          call concoct_ZZ(nt(ipl),fun_mf,nb,np,vmf,dwgh,aa,ds,zz)
c
          ip=ip+np
        enddo
c
        deallocate(dwgh,vmf)
        deallocate(ddif,aa)
c
        return
        end