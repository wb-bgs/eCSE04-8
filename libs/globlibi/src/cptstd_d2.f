cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cptstd_d2
c                V.Lesur 17.08.2011
c
c       Derived from cptstd_2 for diagonal covariance matrix only
c
c       Compute the STD
c
c       input:
c          npmax          number max of data points handled together
c          ipg            number of the starting point
c          npt            number of data points
c          jcov           integer arrays describing cov format
c          cov            Covariance matrix in SLAP column format
c          ddat           data vector
c          xyzf           result of forward modelling
c          fun_std        std function
c
c       output:
c          std            STD value
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cptstd_d2(npmax, ipg, npt, jcov, cov,
     >                       ddat, xyzf, fun_std, std)
c
        implicit none
c
        integer ip,ipg,npt,np,npmax,jcov(*),i,ipl
        real*8 std,ddat(*),xyzf(*),cov(*),std0
        real*8, allocatable :: dwgh(:),ddif(:)
c
        real*8 fun_std
        external fun_std
c
        allocate(dwgh(1:npmax))
        allocate(ddif(1:npmax))
c
        std=0.0d0
c
        ip=1
        do while (ip.le.npt)
          ipl=ipg+ip-1
          np=min0(npt-ip+1,npmax)
c
c  Calculate the inverse covariance matrix
          do i=1,np
            dwgh(i)=1.d0/cov(jcov(ipl+i-1))
          enddo
c
c  calculate the delta data
          do i=1,np
            ddif(i)=ddat(ipl+i-1)-xyzf(ipl+i-1)
          enddo
c
c  Calculate STD
          std0=std
          std=fun_std(ip-1,np,std0,ddif,dwgh)
c
          ip=ip+np
        enddo
c
        deallocate(ddif)
        deallocate(dwgh)
c
        return
        end