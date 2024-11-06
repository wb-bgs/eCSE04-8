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
c       Called: MPI_ALLREDUCE
c
c       input:
c          npts           Total number of points (data + sampling) for all ranks
c          nlocpts        Total number of points for this rank
c          cov            Covariance matrix in SLAP column format
c         j cov           Integer vector describing cov format
c          ddat           data vector
c          xyzf           result of forward modelling
c
c       output:
c          std            STD value
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cptstd_dp(npts, nlocpts, cov,
     >                       jcov, ddat, xyzf, std)
c
        implicit none
c
        include 'mpif.h'
c
        integer npts, nlocpts
        real*8  cov(1:nlocpts)
        integer jcov(1:nlocpts+2)
        real*8  ddat(1:nlocpts)
        real*8  xyzf(1:nlocpts)       
        real*8  std
c
        integer i, ierr
        real*8 dwgh, ddif
c
c
        std = 0.0d0
c
        do i=1,nlocpts
c
c  Calculate the inverse covariance matrix
          dwgh = 1.d0/cov(jcov(i))          
c
c  Calculate the delta data
          ddif = ddat(i)-xyzf(i)

c  Calculate STD
          std = std + dwgh*(ddif**2)
c
        enddo
c
c  Compute the standard deviation across all ranks
        std = dsqrt(std/dble(nlocpts))
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