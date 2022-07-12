ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       subroutine build_damp_space.f   V.Lesur
c
c       Setup equations for the Norm:
c        int_{\Omega} ( B_r \dot B.cm4)^2 d_{\omega}
c
c       input:
c       nub         integer    index for this kind of data
c       npts        integer    current number of data points available
c       nsampts     integer    current number of sampling points
c       nd          integer    such that nd+1 is the lead dim of ddat
c       proc_nsp(*) integer    block lengths
c       proc_isp(*) integer    data pointers
c       ncoeffs     integer    number of coefficients   
c       llm         integer    maxium degree for lithospheric field
c       wgh         real*8     weight for damping
c       bc          real*8     array of initial parameters
c
c       output:
c       nt(*)       integer    array seting data type to nub
c       ddat(*)     real*8     array of data values
c       ijcov(*)    integer    array defining cov matrix
c       cov(*)      real*8     array of cov matrix
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine build_damp_space(nub, npts, nsampts, nd,
     >                              proc_nsp, proc_isp, 
     >                              ncoeffs, llm, wgh, bc,
     >                              nt, ijcov, cov, ddat)
c
        use sph_wmam
c
        implicit none
c
        include 'mpif.h'
c
        integer nub, npts, nsampts
        integer nd, ncoeffs, llm, ijcov(npts+nsampts,*), nt(*)
        integer proc_nsp(*), proc_isp(*)
        real*8 bc(*), ddat(nd+1,*), cov(*), wgh
c
        real*8, allocatable :: vrt(:,:)
        real*8, allocatable :: glw(:)
        real*8, allocatable :: dw(:)
        real*8, allocatable :: dw2(:,:)
c
        integer i,j
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
c
c  Define sampling points
        allocate (glw(nsampts))
        allocate (vrt(nd+1,nsampts))
        allocate (dw2(2,nsampts))
        call set_FG_sampling(llm,nsampts,dw2,glw)
        do i=1,nsampts
          vrt(1,i)=dw2(1,i)
          vrt(2,i)=dw2(2,i)
          vrt(3,i)=6371.2d0
          vrt(4,i)=ryg
        enddo
        deallocate (dw2)
c
        allocate(dw(1:nsampts))
c
        nt(npts+1:npts+nsampts)=1
        dw=0.0d0
        if (rank.eq.0) write(*,*) ' Calculating CM4 components'
        call cpt_dat_vals_p(nd,proc_nsp,proc_isp,nt(npts+1),
     >                      vrt,ncoeffs,bc,sph_bi,dw)
        vrt(5,1:nsampts)=dw(1:nsampts)
        if (rank.eq.0) write(*,*) '  X CM4 component calculated'
c
        nt(npts+1:npts+nsampts)=2
        dw=0.0d0
        call cpt_dat_vals_p(nd,proc_nsp,proc_isp,nt(npts+1),
     >                      vrt,ncoeffs,bc,sph_bi,dw)
        vrt(6,1:nsampts)=dw(1:nsampts)
        if (rank.eq.0) write(*,*) '  Y CM4 component calculated'
c
        nt(npts+1:npts+nsampts)=3
        dw=0.0d0
        call cpt_dat_vals_p(nd,proc_nsp,proc_isp,nt(npts+1),
     >                      vrt,ncoeffs,bc,sph_bi,dw)
        vrt(7,1:nsampts)=dw(1:nsampts)
        if (rank.eq.0) then
          write(*,*) '  Z CM4 component calculated'
          write(*,*) ''
        endif

        deallocate(dw)
c
c for each data point define the covariance
        j=npts
        do i=1,nsampts
          j=j+1
          ddat(1:nd+1,j)=0.0d0
          ddat(1:7,j)=vrt(1:7,i)
          nt(j)=nub
c
          ijcov(j,1)=j
          ijcov(j,2)=j
          cov(j)=1.d0/glw(i)
          cov(j)=cov(j)/wgh
        enddo
c
        deallocate (vrt)
        deallocate (glw)
c
        return
        end subroutine build_damp_space