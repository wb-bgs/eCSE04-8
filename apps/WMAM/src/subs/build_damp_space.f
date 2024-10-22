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
c       ryg             real*8     reference year for the model
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
     >                              nd, ncoeffs, llm, ryg, wgh, bc,
     >                              ijcov, cov, ppos)
c
        implicit none
c
        include 'mpif.h'
c
        real*8, parameter :: RAG = 6371.2d0
c
        integer nlocdatpts, nlocsampts
        integer imin_locpts, imin_locsampts
        integer nd, ncoeffs, llm, ijcov(nlocdatpts+nlocsampts+2,*)
        real*8 bc(*), ppos(nd+1,*), cov(*), ryg, wgh
c
        real*8, allocatable :: vrt(:,:)
        real*8, allocatable :: glw(:)
        real*8, allocatable :: cm(:)
        real*8, allocatable :: dw(:,:)
c
        integer i, j, k
        integer ierr,rank
c
c        
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
c  Define sampling points
        allocate(glw(nlocsampts))
        allocate(vrt(nd+1,nlocsampts))
        allocate(dw(2,nlocsampts))
        call set_FG_sampling(llm, imin_locsampts, 
     >                       imin_locsampts+nlocsampts-1,
     >                       dw, glw)
        do i = 1,nlocsampts
          vrt(1,i) = dw(1,i)
          vrt(2,i) = dw(2,i)
          vrt(3,i) = RAG
          vrt(4,i) = ryg
        enddo
        deallocate(dw)
c
c  Calculate CM4 components 
        if (rank .eq. 0) write(*,*) 'Calculating CM4 components'
        allocate(cm(ncoeffs))
        do i = 1,3
          j = 4+i
          cm = 0.0d0 
          do k = 1,nlocsampts
            call sph_bi('f', i, ncoeffs, bc, vrt(1,k), cm)
            vrt(j,k) = cm(1)
          enddo
        enddo
        deallocate(cm)
        if (rank .eq. 0) write(*,*) 'XYZ CM4 components calculated'

c
c  For each data point define the covariance
        j = nlocdatpts+1
        k = imin_locpts+nlocdatpts
        do i = 1,nlocsampts
          ppos(1:nd,j) = vrt(1:nd,i)
          ppos(nd+1,j) = 0.0d0
c
          ijcov(j,1) = k
          ijcov(j,2) = k
          cov(j) = 1.d0/glw(i)
          cov(j) = cov(j)/wgh

          j = j+1
          k = k+1
        enddo
c
        deallocate(vrt)
        deallocate(glw)
c
        return
        end subroutine build_damp_space