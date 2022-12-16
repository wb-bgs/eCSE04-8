ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       subroutine build_damp_space.f   V.Lesur
c
c       Setup equations for the Norm:
c        int_{\Omega} ( B_r \dot B.cm4)^2 d_{\omega}
c
c       input:
c       nlocdatpts      integer    number of data points local to rank
c       nlocsampts      integer    number of sampling points local to rank
c       imin_locpts     integer    local rank's global index for data+sampling points
c       imin_locsampts  integer    local rank's global index for sampling points
c       nd              integer    such that nd+1 is the lead dim of ppos
c       ncoeffs         integer    number of coefficients   
c       llm             integer    maxium degree for lithospheric field
c       wgh             real*8     weight for damping
c       bc              real*8     array of initial parameters
c
c       output:
c       ijcov(*)        integer    array defining cov matrix
c       cov(*)          real*8     array of cov matrix
c       ppos(*)         real*8     array of data values
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine build_damp_space(nlocdatpts, nlocsampts,
     >                              imin_locpts, imin_locsampts,
     >                              nd, ncoeffs, llm, wgh, bc,
     >                              ijcov, cov, ppos)
c
        use sph_wmam
c
        implicit none
c
        include 'mpif.h'
c
        integer nlocdatpts, nlocsampts
        integer imin_locpts, imin_locsampts
        integer nd, ncoeffs, llm, ijcov(nlocdatpts+nlocsampts+2,*)
        real*8 bc(*), ppos(nd+1,*), cov(*), wgh
c
        real*8, allocatable :: vrt(:,:)
        real*8, allocatable :: glw(:)
        real*8, allocatable :: dw2(:)
        real*8, allocatable :: dw3(:,:)
c
        integer i,j,k
        real*8 r2,qpi
c
        external sph_bi
c
c  MPI, determine rank
        integer ierr,rank
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
        r2=rag**2
        qpi=16.d0*datan(1.d0)
c
c  Define sampling points
        allocate (glw(nlocsampts))
        allocate (vrt(nd+1,nlocsampts))
        allocate (dw3(2,nlocsampts))
        call set_FG_sampling(llm, imin_locsampts, 
     >                       imin_locsampts+nlocsampts-1,
     >                       dw3,glw)
        do i=1,nlocsampts
          vrt(1,i)=dw3(1,i)
          vrt(2,i)=dw3(2,i)
          vrt(3,i)=rag
          vrt(4,i)=ryg
        enddo
        deallocate (dw3)
c
        allocate(dw2(1:nlocsampts))
c
        dw2=0.0d0
        if (rank.eq.0) write(*,*) ' Calculating CM4 components'
        call cpt_dat_vals_p(nd, nlocsampts, 1, vrt, ncoeffs,
     >                      bc, sph_bi, dw2)
        vrt(5,1:nlocsampts)=dw2(1:nlocsampts)
        if (rank.eq.0) write(*,*) '  X CM4 component calculated'
c
        dw2=0.0d0
        call cpt_dat_vals_p(nd, nlocsampts, 2, vrt, ncoeffs,
     >                      bc, sph_bi, dw2)
        vrt(6,1:nlocsampts)=dw2(1:nlocsampts)
        if (rank.eq.0) write(*,*) '  Y CM4 component calculated'
c
        dw2=0.0d0
        call cpt_dat_vals_p(nd, nlocsampts, 3, vrt, ncoeffs,
     >                      bc, sph_bi, dw2)
        vrt(7,1:nlocsampts)=dw2(1:nlocsampts)
        if (rank.eq.0) then
          write(*,*) '  Z CM4 component calculated'
          write(*,*) ''
        endif

        deallocate(dw2)
c
c for each data point define the covariance
        j=nlocdatpts+1
        k=imin_locpts+nlocdatpts
        do i=1,nlocsampts
          ppos(1:nd,j)=vrt(1:nd,i)
          ppos(nd+1,j)=0.0d0
c
          ijcov(j,1)=k
          ijcov(j,2)=k
          cov(j)=1.d0/glw(i)
          cov(j)=cov(j)/wgh

          j=j+1
          k=k+1
        enddo
c
        deallocate (vrt)
        deallocate (glw)
c
        return
        end subroutine build_damp_space
