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
c          npm            Maximum total number of data points
c          npt            total number of data points
c          ipg            where to start in data file!
c          nd             space dimension
c          ppos           data point position in ndD
c          nb             Number or base function to use
c          FM             misfit function (like l2_norm.f)
c          BS             the "Base functions" subroutine to use
c          bc             Estimation of Base function coefficients
c          ds             current descent direction
c          icov/jcov      integer arrays describing cov format
c          cov            Covariance matrix in SLAP column format
c          ddat           data vector
c          nt             vector indicating data type
c          xyzf           result of forward modelling
c
c       output:
c          ZZ             Vector A^t.W.A.ds (nb)
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cptAtWAds_p(npmax,npm,npt,ipg,nd,ppos,nb,FM,BS,bc,ds
     >                        ,icov,jcov,cov,ddat,nt,xyzf,ZZ)
c
        implicit none
c
        include 'mpif.h'
c
        integer npt,npmax,npm,ipg,nd,nb,icov(*),jcov(*),nt(*)
        real*8 ddat(*),xyzf(*),cov(*),ppos(*),bc(*),ds(*)
        real*8 ZZ(*)
c
        integer :: i,ip,np,inp
        real*8, allocatable :: ZZt(:)
c
        real*8 FM
        external FM,BS
c
c  All defining parallel enviroment
        integer ierr,rank,size
        integer, allocatable :: proc_np(:),proc_ip(:)
        call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
        allocate(proc_np(1:size),proc_ip(1:size))
c
c  All: define what is going to each SP
        i=npt-ipg+1
        call thread_segmenter(size,i,proc_np,proc_ip)
        proc_ip(1:size)=proc_ip(1:size)+ipg-1
c
c  All: find out what are the data to work with
        np=proc_np(rank+1)
        ip=proc_ip(rank+1)
c
c  All: Now does the work
        allocate(ZZt(1:nb))
        call cptAtWAds(npmax,npm,np,ip,nd,ppos,nb,FM,BS,bc,ds
     >                   ,icov,jcov,cov,ddat,nt,xyzf,ZZt)
c
c  All: Gather & SUM the ZZt results from ALL the other Processes
        call MPI_ALLREDUCE(ZZt,ZZ,nb,MPI_DOUBLE_PRECISION
     >                            ,MPI_SUM,MPI_COMM_WORLD,ierr)
c
        deallocate(ZZt)
        deallocate(proc_np,proc_ip)
        return
        end
