ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       subroutine build_damp_space.f   V.Lesur
c
c       Setup equations for the Norm:
c        int_{\Omega} ( B_r \dot B.cm4)^2 d_{\omega}
c
c       input:
c       nub       integer    index for this kind of data
c       npm       integer    maximum number of data points acceptable
c       npt       integer    current number of data points available
c       ncm       integer    maximum number of covarience values
c       nc        integer    current number of covarience values
c       nd        integer    such that nd+1 is the lead dim of ddat
c       llm       integer    maxium degree for lithospheric field
c       wgh       real*8     weight for damping
c
c       output:
c       npt       integer    current number of data points available
c       nc        integer    current number of covariance values
c       nt(*)     integer    array seting data type to nub
c       ddat(*)   real*8     array of data values
c       ijcov(*)  integer    array defining cov matrix
c       cov(*)    real*8     array of cov matrix
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine build_damp_space(nub,npm,npt,ncm,nc,nd
     >                                   ,llm,wgh,nt,ijcov,cov,ddat)
c
        use sph_wmam
c
        implicit none
c
        include 'mpif.h'
c
        integer nub,npm,npt,ncm,nc,nd,llm,ijcov(ncm,*),nt(*)
        real*8 ddat(nd+1,*),cov(*),wgh
c
        integer, parameter :: nbm=15*17 
        integer nb(2)
        real*8 bc(nbm)
        character nom*100
c
        integer np
        real*8, allocatable :: vrt(:,:)
        real*8, allocatable :: glw(:)
        real*8, allocatable :: dw(:)
        real*8, allocatable :: dw2(:,:)
c
        integer i
        real*8 r2,qpi
c
        external sph_bi
c
c  MPI, determine rank
        integer ierr,rank
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
        r2=6371.2d0
        r2=r2**2
        qpi=16.d0*datan(1.d0)
        if(rank.eq.0)write(*,*)'4 Pi = ',qpi
c
c  Define sampling points
        np=(llm+1)*(2*llm+1)
        allocate (glw(np))
        allocate (vrt(nd+1,np))
        allocate (dw2(2,np))
        call set_FG_sampling(llm,np,dw2,glw)
        do i=1,np
          vrt(1,i)=dw2(1,i)
          vrt(2,i)=dw2(2,i)
          vrt(3,i)=6371.2d0
          vrt(4,i)=ryg
        enddo
        deallocate (dw2)
c
c  Define the reference field at these points
        allocate(dw(1:np))
        if(rank.eq.0)write(*,*)"define reference field values"
        nom='./Data/coef_1990_15.dat'
c
        nb(1)=nbm
        i=nbm
        call read_model(nom,ryg,bc,i)
	if(rank.eq.0)write(*,*)"have read file"
c
        if(rank.eq.0)write(*,*) npt, np
        nt(npt+1:npt+np)=1
        dw=0.0d0
	if(rank.eq.0)write(*,*)"entering cpt_dat_vals_p"
        call cpt_dat_vals_p(nd,np,np,nt(npt+1),vrt,nb(1),bc,sph_bi,dw)
        vrt(5,1:np)=dw(1:np)
        if(rank.eq.0)write(*,'(A)')"X CM4 component calculated"
c
        nt(npt+1:npt+np)=2
        dw=0.0d0
        call cpt_dat_vals_p(nd,np,np,nt(npt+1),vrt,nb(1),bc,sph_bi,dw)
        vrt(6,1:np)=dw(1:np)
        if(rank.eq.0)write(*,'(A)')"Y CM4 component calculated"
c
        nt(npt+1:npt+np)=3
        dw=0.0d0
        call cpt_dat_vals_p(nd,np,np,nt(npt+1),vrt,nb(1),bc,sph_bi,dw)
        vrt(7,1:np)=dw(1:np)
        if(rank.eq.0)write(*,'(A)')"Z CM4 component calculated"
        deallocate(dw)
c
c for each data point define the covariance
        do i=1,np
          npt=npt+1
          ddat(1:nd+1,npt)=0.0d0
          ddat(1:7,npt)=vrt(1:7,i)
          nt(npt)=nub
c
          nc=nc+1
          ijcov(nc,1)=npt
          ijcov(nc,2)=npt
          cov(nc)=1.d0/glw(i)
c         cov(nc)=cov(nc)*r2*qpi
          cov(nc)=cov(nc)/wgh
        enddo
c
        deallocate (vrt)
        deallocate (glw)
c
        return
        end subroutine build_damp_space
