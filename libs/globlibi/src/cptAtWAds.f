cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cptAtWAds
c		Vincent Lesur  18/08/2011
c
c      derived from ssqgh_d.f, for diagonal COV matrix  only
c       Should be used with F90 ot later versions
c
c       input:
c          npmax          number max of data points handled together
c          nlocpts        number of data+sampling points local to rank
c          nlocdatpts    number of data points assigned to rank
c          shdeg         max SH degree value
c          ipg            where to start in data file!
c          nd             space dimension
c          ppos           data point position in ndD
c          nb             Number or base function to use
c          bc             Estimation of Base function coefficients
c          ds             current descent direction
c          jcov           integer arrays describing cov format
c          cov            Covariance matrix in SLAP column format
c          ddat           data vector
c          xyzf           result of forward modelling
c
c       output:
c          zz             Vector A^t.W.A.ds (nb)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cptAtWAds(npmax, nlocpts, nlocdatpts, shdeg, ipg,
     >                       nd, d2a, ppos, nb,
     >                       bc, ds,
     >                       jcov, cov, ddat,
     >                       xyzf, zz)
c
        implicit none
c
        integer ip,nlocpts,np,npmax,nd,nb,jcov(*)
        integer i,ipg,ipl,nlocdatpts,shdeg
        real*8 d2a(*),ddat(*),xyzf(*),cov(*),ppos(nd+1,*),bc(*),ds(*)
        real*8 zz(*)
c       real*8, allocatable :: vmf(:)
        real*8, allocatable :: dwgh(:),ddif(:),aa(:,:)
c
c
c  ipg : ip global
c  ipl : ip local
c
c       allocate(vmf(1:npmax))
        allocate(dwgh(1:npmax),ddif(1:npmax))
        allocate(aa(nb,npmax))
c
        zz(1:nb)=0.0d0
c
        ip=1
        do while (ip.le.nlocpts) 
          ipl=ipg+ip-1
          np=min0(nlocpts-ip+1,npmax)
c
c  calculate the equations of condition
          call mkArows(shdeg, nb, nd, np,
     >                 ip, nlocdatpts,
     >                 d2a, bc, ppos(1,ip),
     >                 aa)

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
c           vmf(i)=ddif(i)*dsqrt(dwgh(i))
c         enddo
c
c  update the G matrix and B vector
          call concoct_ZZ(nb,np,dwgh,aa,ds,zz)
c
          ip=ip+np
        enddo
c
c       deallocate(vmf)
        deallocate(dwgh,ddif,aa)
c
        return
        end