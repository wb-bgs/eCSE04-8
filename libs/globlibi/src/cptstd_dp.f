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
c modified 14.06.20011 (V. Lesur)
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
c          npt            total number of data points
c          nt             array of data type
c          icov/jcov      integer arrays describing cov format
c          cov            Covariance matrix in SLAP column format
c          ddat           data vector
c          xyzf           result of forward modelling
c          fstd           std function
c
c       output:
c          std            STD value
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cptstd_dp(npmax,npt,nt,icov,jcov,cov
     >                                            ,ddat,xyzf,fstd,std)
c
        implicit none
c
        include 'mpif.h'
c
        integer :: npmax,npt,icov(*),jcov(*),nt(*)
        integer :: i,ip,np,inp,irank,status(MPI_STATUS_SIZE)
        integer, allocatable :: proc_np(:),proc_ip(:)
        real*8  :: cov(*),ddat(*),xyzf(*),std,dw
        real*8, allocatable :: vstd(:),vnp(:)
c
        real*8 fstd
        external fstd
c
        integer ierr,rank,size
        call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
        allocate(proc_np(1:size),proc_ip(1:size))
        allocate(vstd(1:size),vnp(1:size))
c
c  All: define what is going to each SP 
        i=npt
        call thread_segmenter(size,i,proc_np,proc_ip)
c
c  All: find out what are the data to work with
        np=proc_np(rank+1)
        ip=proc_ip(rank+1)
c
c  All: Now does the work
        call cptstd_d2(npmax,ip,np,nt,icov,jcov,cov,ddat,xyzf
     >                ,fstd,std)
c
c  All: Gather the results from other Processes
        call MPI_ALLGATHER(std,1,MPI_DOUBLE_PRECISION,vstd
     >       ,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
c
c  ALL: put together all STDs
        vnp=dble(proc_np)
        np=-size
        ip=0
        std=fstd(nt,ip,np,dw,vstd,vnp)
c
        deallocate(vstd,vnp)
        deallocate(proc_np,proc_ip)
        return
        end
