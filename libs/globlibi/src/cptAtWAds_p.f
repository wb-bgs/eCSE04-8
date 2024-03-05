cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cptAtWAds_p
c		Vincent Lesur 19/08/2011
c
c       Parallel interface for cptAtWAds.f
c
c       Called: cptAtWAds, MPI_ALLREDUCE
c
c       input:
c          npmax          number max of data point with correlated errors
c          nd             space dimension
c          nlocpts        number of data+sampling points local to rank
c          nlocdatpts    number of data points assigned to rank
c          shdeg         max SH degree value
c          d2a           pre-computed array for mklf_F2()
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
        subroutine cptAtWAds_p(npmax, nd, nlocpts, nlocdatpts, shdeg,
     >                         d2a, ppos, nb, bc, ds,
     >                         jcov, cov, ddat,
     >                         xyzf, zz)
c
        implicit none
c
        include 'mpif.h'
c
        integer npmax,nd,nb,jcov(*)
        integer nlocpts,nlocdatpts,shdeg
        real*8 d2a(*),ddat(*),xyzf(*),cov(*),ppos(*),bc(*),ds(*)
        real*8 zz(*)
c
        integer :: ip,np
c
c
c  All defining parallel enviroment
        integer ierr,rank
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
c  All: Now does the work
        call cptAtWAds(npmax, nlocpts, nlocdatpts, shdeg, 1,
     >                 nd, d2a, ppos, nb,
     >                 bc, ds,
     >                 jcov, cov, ddat,
     >                 xyzf, zz)

c
c  All: Gather & SUM the zzt results from ALL the other Processes
        call MPI_ALLREDUCE(MPI_IN_PLACE, zz, nb, MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
        return
        end