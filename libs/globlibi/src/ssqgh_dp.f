cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine ssqgh_dp
c		Vincent Lesur 13/06/2006
c
c       Modified 18.08.2011 to call ssqgh_d in place of ssqgh as
c       none diagonal covariance matrix element are not expected
c       for Large system (V.Lesur) 
c
c       Parallel interface for ssqgh.f
c
c       Called: ssqgh, MPI_ALLREDUCE
c
c       input:
c          npmax          number max of data point with correlated errors
c          nd             space dimension
c          nlocdatpts     number of data points local to rank
c          nlocpts        number of data+sampling points local to rank
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
        subroutine ssqgh_dp(npmax, nd, nlocdatpts, nlocpts,
     >                      ppos, nb, fun_mf, sub_base, bc,
     >                      jcov, cov, ddat,
     >                      xyzf, gj, hj)
c
        implicit none
c
        include 'mpif.h'
c
        integer npmax,nd,nb,jcov(*)
        integer nlocdatpts,nlocpts
        real*8 ddat(*),xyzf(*),cov(*),ppos(*),bc(*)
        real*8 gj(*),hj(*)
c
        integer :: nb2
        real*8, allocatable :: gj_hj(:)
c
        real*8 fun_mf
        external fun_mf,sub_base
c
c  All defining parallel enviroment
        integer ierr,rank
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
c  All: Now does the work
        nb2=nb*2
        allocate(gj_hj(nb2))
c
        call ssqgh_d(npmax, nlocdatpts, nlocpts, 1,
     >               nd, ppos, nb,
     >               fun_mf, sub_base, bc,
     >               jcov, cov,
     >               ddat, xyzf,
     >               gj_hj(1), gj_hj(nb+1))
c
c  All: Gather & SUM the GJ and HJ results from ALL the other Processes
        call MPI_ALLREDUCE(MPI_IN_PLACE, gj_hj, nb2,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
        gj(1:nb)=gj_hj(1:nb)
        hj(1:nb)=gj_hj(nb+1:nb2)
c
        deallocate(gj_hj)
c
        return
        end