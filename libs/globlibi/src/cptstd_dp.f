cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cptstd_dp
c		Vincent Lesur 28/04/2005
c
c modified for diagonal cov matrix only (lesur 18.08.2011)
c
c modified 14.07.2011 (V. Lesur)
c          2- data type nt(*) included as input parameter
c          1- call to cptstd_2 updated
c
c modified 14.06.2011 (V. Lesur)
c       set ip=proc_ip(rank+1) in place of proc_ip(rank+1)+1
c       (bug introduced at last modification)
c
c modified for subroutine thread_segmenter.f (lesur 11.06.2010)
c
c modified 12.09.2007 V.Lesur
c       replace common by MPI_Comm_size and MPI_Comm_rank
c
c       Parallel interface for cptstd.f
c
c       Called: cptstd_d2, MPI_ALLGATHER
c
c       input:
c          npmax          number max of data point with correlated errors
c          npts           Total number of points (data + sampling) for all ranks
c          nlocpts        Total number of points for this rank
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
        subroutine cptstd_dp(npmax, npts, nlocpts, jcov, cov, ddat,
     >                       xyzf, fun_std, std)
c
        implicit none
c
        include 'mpif.h'
c
        integer :: npmax,npts,nlocpts,jcov(*)
        real*8  :: cov(*),ddat(*),xyzf(*)       
        real*8 fun_std
        external fun_std
        real*8 std

        integer ierr
c
c
c  All: Now does the work
        call cptstd_d2(npmax, 1, nlocpts, jcov,
     >                 cov, ddat, xyzf,
     >                 fun_std, std)
c
c  Compute the standard deviation across all ranks
        std = nlocpts*(std**2)
c
        call MPI_ALLREDUCE(MPI_IN_PLACE, std, 1,
     >                     MPI_DOUBLE_PRECISION,
     >                     MPI_SUM, MPI_COMM_WORLD, ierr)
c
        std = dsqrt(std/dble(npts))
c        
        return
        end