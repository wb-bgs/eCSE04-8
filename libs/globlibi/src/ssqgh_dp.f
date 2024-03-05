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
c          nlocpts        number of data+sampling points local to rank
c          nlocdatpts    number of data points assigned to rank
c          shdeg         max SH degree value
c          d2a           pre-computed array for mklf_F2()
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
        subroutine ssqgh_dp(npmax, nd, nlocpts, nlocdatpts, shdeg,
     >                      d2a, ppos, nb, bc,
     >                      jcov, cov, ddat,
     >                      xyzf, gj, hj)
c
        implicit none
c
        include 'mpif.h'
c
        integer npmax,nd,nb,jcov(*)
        integer nlocpts,nlocdatpts,shdeg
        real*8 d2a(*),ddat(*),xyzf(*),cov(*),ppos(*),bc(*)
        real*8 gj(*),hj(*)
c
c
c  All defining parallel enviroment
        integer ierr,rank
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
c  All: Now does the work
        call ssqgh_d(npmax, nlocpts, nlocdatpts, shdeg,
     >               d2a, nd, ppos, nb, bc,
     >               jcov, cov,
     >               ddat, xyzf,
     >               gj, hj)
c
c  All: Gather & SUM the GJ and HJ results from ALL the other Processes
        call MPI_ALLREDUCE(MPI_IN_PLACE, gj, nb,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
        call MPI_ALLREDUCE(MPI_IN_PLACE, hj, nb,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
        return
        end