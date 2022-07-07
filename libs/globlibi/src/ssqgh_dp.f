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
c          npm            Maximum total number of data points
c          npt            total number of data points
c          ipg            where to start in data file!
c          nd             space dimension
c          ppos           data point position in ndD
c          nb             Number or base function to use
c          FM             misfit function (like l2_norm.f)
c          BS             the "Base functions" subroutine to use
c          bc             Estimation of Base function coefficients
c          icov/jcov      integer arrays describing cov format
c          cov            Covariance matrix in SLAP column format
c          ddat           data vector
c          nt             vector indicating data type
c          xyzf           result of forward modelling
c
c       output:
c          GJ             gradient of the weighted sum of squares (nb)
c          Hj             diagonal of the Hessian (nb)
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine ssqgh_dp(npmax,npm,npt,ipg,nd,ppos,nb,FM,BS,bc
     >                        ,icov,jcov,cov,ddat,nt,xyzf,GJ,HJ)
c
        implicit none
c
        include 'mpif.h'
c
        integer npt,npmax,npm,ipg,nd,nb,icov(*),jcov(*),nt(*)
        real*8 ddat(*),xyzf(*),cov(*),ppos(*),bc(*)
        real*8 GJ(*),HJ(*)
c
        integer :: i,ip,np
        real*8, allocatable :: GJt(:),HJt(:)
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
        allocate(GJt(1:nb),HJt(1:nb))
        call ssqgh_d(npmax,npm,np,ip,nd,ppos,nb,FM,BS,bc
     >                   ,icov,jcov,cov,ddat,nt,xyzf,GJt,HJt)
c
c  All: RE-scale GJt and HJt for npt points
c       GJt(1:nb)=dble(np)*GJt(1:nb)/dble(npt)
c       HJt(1:nb)=dble(np)*HJt(1:nb)/dble(npt)
c
c  All: Gather & SUM the GJ results from ALL the other Processes
        call MPI_ALLREDUCE(GJT,GJ,nb,MPI_DOUBLE_PRECISION
     >                            ,MPI_SUM,MPI_COMM_WORLD,ierr)
c
c  All: Make sure that it's all & well done (may not be necessary)
c       call MI_BARRIER(MPI_COMM_WORLD,ierr)
c
c  All: Gather & SUM the HJ results from ALL the other Processes
        call MPI_ALLREDUCE(HJT,HJ,nb,MPI_DOUBLE_PRECISION
     >                            ,MPI_SUM,MPI_COMM_WORLD,ierr)
c
c  All: Make sure that it's all & well done (may not be necessary)
c       call MI_BARRIER(MPI_COMM_WORLD,ierr)
c
        deallocate(GJt,HJt)
        deallocate(proc_np,proc_ip)
        return
        end
