cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Program MOD_WMAM
c    
c    Build a spherical harmonic model from a regular grid of 
c    total intensity data in Geocentric
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	program mod_wmam
c
        use sph_wmam
c
        implicit none
c
        include 'mpif.h'
c
        character(len=*), parameter :: VERSION="4.2.0"
c
        integer, parameter :: POLAK_RIBIERE=1
        integer, parameter :: CONJUGATE_GRADIENT=2
        integer, parameter :: ND=7
        integer, parameter :: NPMAX=3
c
        character fname*100
        character buf*100
        integer fhandle
        integer i, j, ix, iy, il, im
        integer nx, ny
        integer itmax(3), nub
        real*8  std, stdt, dd, dl(3), dampfac
        real*8  resdeg
        integer cmdcnt, shdeg, scheme, serialrd
        integer ncoeffs, nparams
        integer ndatpts, nsampts, npts
        integer nlocdatpts, imin_locdatpts
        integer nlocsampts, imin_locsampts
        integer nlocpts, imin_locpts
c
        character(100) :: argstr
c
        integer, allocatable :: proc_ndp(:), proc_idp(:)
        integer, allocatable :: proc_nsp(:), proc_isp(:)
        integer, allocatable :: proc_np(:), proc_ip(:)
c
        integer, allocatable :: ijcov(:,:)
        real*8, allocatable :: ppos(:,:), bc(:), cov(:), dw(:)
        real*8, allocatable :: gg(:,:), bb(:)
        real*8, allocatable :: err(:)
        real*8 diff(4)
c
c  globlibi functions
        real*8 dsind, l2_std, l2_norm
        external dsind, l2_std, l2_norm
c
c  globlibi subroutines
        external sph_bi, damp_rien

c
c  Initialize MPI, determine rank
        integer ierr, nranks, rank
        call MPI_Init(ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, nranks, ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
c  Read in command line arguments
        cmdcnt = COMMAND_ARGUMENT_COUNT()
        if (cmdcnt.ge.1) then
          call GET_COMMAND_ARGUMENT(1,argstr)
          read(argstr,*) shdeg
        else
          shdeg=200
        endif
        if (cmdcnt.ge.2) then
          call GET_COMMAND_ARGUMENT(2,argstr)
          read(argstr,*) resdeg
        else
          resdeg=1.0
        endif
        if (cmdcnt.ge.3) then
          call GET_COMMAND_ARGUMENT(3,argstr)
          read(argstr,*) scheme
        else
          scheme=POLAK_RIBIERE
        endif
        if (cmdcnt.ge.4) then
          call GET_COMMAND_ARGUMENT(4,argstr)
          read(argstr,*) dampfac
        else
          dampfac=5.0d0
        endif
        if (cmdcnt.ge.5) then
          call GET_COMMAND_ARGUMENT(5,argstr)
          read(argstr,*) serialrd
          if (serialrd .lt. 0) then
            serialrd = 0
          elseif (serialrd .gt. 2) then
            serialrd = 2
          endif
        else
          serialrd=0
        endif
        if (rank.eq.0) then
          write(*,*) 'WMAM v', VERSION
          write(*,*) ''
          write(*,*) 'shdeg: ', shdeg
          write(*,*) 'resdeg: ', resdeg
          write(*,*) 'scheme: ', scheme
          write(*,*) 'dampfac: ', dampfac
          write(*,*) 'serialrd: ', serialrd
          write(*,*) ''
        endif

c
c  Settings
        nparams=shdeg*(shdeg+2)
        nx=nint(1.0/resdeg)*180-1
        ny=nint(1.0/resdeg)*360
        ndatpts=nx*ny
        nsampts=(shdeg+1)*(2*shdeg+1)
        ncoeffs=255
        npts=ndatpts+nsampts

        if (rank.eq.0) then
          write(*,*) 'MPI Ranks:', nranks
          write(*,*) 'Parameters: ', nparams
          write(*,*) 'nx: ', nx
          write(*,*) 'ny: ', ny
          write(*,*) 'Data points: ', ndatpts
          write(*,*) 'Sampling points: ', nsampts
          write(*,*) 'Data+Sampling points: ', npts
          write(*,*) ''
        endif

c
c  Partition workload
        allocate(proc_ndp(nranks),proc_idp(nranks))
        call thread_segmenter(nranks,ndatpts,proc_ndp,proc_idp)
        nlocdatpts = proc_ndp(rank+1)
        imin_locdatpts = proc_idp(rank+1)
        
        allocate(proc_nsp(nranks),proc_isp(nranks))        
        call thread_segmenter(nranks,nsampts,proc_nsp,proc_isp)
        nlocsampts = proc_nsp(rank+1)
        imin_locsampts = proc_isp(rank+1)

        allocate(proc_np(nranks),proc_ip(nranks))
c  initialise proc_np and proc_ip arrays, so that the vnp array
c  used by the cptstd_dp subroutine (from the globlibi library)
c  can be setup correctly
        proc_np(1:nranks)=proc_ndp(1:nranks)+proc_nsp(1:nranks)
        proc_ip(1)=1
        do i=2,nranks
          proc_ip(i)=proc_ip(i-1)+proc_np(i-1)
        enddo
        nlocpts = proc_np(rank+1)
        imin_locpts = proc_ip(rank+1)

        deallocate(proc_ndp,proc_idp)
        deallocate(proc_nsp,proc_isp)

c
c  Array allocations
        allocate(bc(nparams))
        allocate(ppos(ND+1,nlocpts))
        allocate(cov(nlocpts))
        allocate(ijcov(nlocpts+2,2))
        allocate(dw(nlocpts))

c
c  Initialize sph_wmam module
        call init_sph_wmam(shdeg, nparams, nlocdatpts)

c
c  Output array sizes
        if (rank.eq.0) then
          write(*,*) ''
          write(*,*) 'MPI Rank:', rank
          write(*,*) ''
          write(*,*) 'Local Data points: ', nlocdatpts
          write(*,*) 'Global Index for Data points: ',
     >               imin_locdatpts
          write(*,*) ''
          write(*,*) 'Local Sampling points: ', nlocsampts
          write(*,*) 'Global Index for Sampling points: ',
     >               imin_locsampts
          write(*,*) ''
          write(*,*) 'Local Data+Sampling points: ', nlocpts
          write(*,*) 'Global Index for Data+Sampling points: ',
     >               imin_locpts
          write(*,*) ''
          write(*,*) ''
        endif

c
c  Read in reference model
        bc(1:nparams)=1.0d0
        fname='./Data/coef_1990_15.dat.bin'
        if (rank.eq.0) write(*,*)
     >    'Reading in reference model, ', fname
        if (serialrd .gt. 0) then
          call mpi_read_ref_model(fname, ncoeffs, bc)
        else
          call mpi_read_all_ref_model(fname, ncoeffs, bc)
        endif
        if (rank.eq.0) then
          write(*,*) 'Coefficients: ', ncoeffs
          write(*,*) ''
        endif
        if (ncoeffs .eq. 0) stop

c
c  Read in data
        fname='./Data/wdmam_geocentric.dat.bin'
        if (rank.eq.0) write(*,*)
     >    'Reading in data, ', fname
        call mpi_read_all_data(fname, ND, ndatpts, nlocdatpts,
     >                         imin_locdatpts, ppos)
        if (nlocdatpts .eq. 0) stop
        
c
c  Calculate CM4 components 
        if (rank.eq.0) write(*,*) 'Calculating CM4 components'

        dw=0.0d0
        call cpt_dat_vals_p(ND, nlocdatpts, 1, ppos, ncoeffs,
     >                      bc, sph_bi, dw)
        ppos(5,1:nlocdatpts)=dw(1:nlocdatpts)
        if (rank.eq.0) write(*,*) ' X CM4 component calculated'
        
        dw=0.0d0
        call cpt_dat_vals_p(ND, nlocdatpts, 2, ppos, ncoeffs,
     >                      bc, sph_bi, dw)
        ppos(6,1:nlocdatpts)=dw(1:nlocdatpts)
        if (rank.eq.0) write(*,*) ' Y CM4 component calculated'

        dw=0.0d0
        call cpt_dat_vals_p(ND, nlocdatpts, 3, ppos, ncoeffs,
     >                      bc, sph_bi, dw)
        ppos(7,1:nlocdatpts)=dw(1:nlocdatpts)
        if (rank.eq.0) then
          write(*,*) ' Z CM4 component calculated'
          write(*,*) ''
        endif

c
c  Define covariance matrix: sin(colat) weight
        if (rank.eq.0) then
          write(*,*) 'Define covariance matrix'
          write(*,*) ''
        endif
        j=imin_locpts
        do i=1,nlocdatpts
          cov(i)=dsind(ppos(1,i))
          cov(i)=1.d0/cov(i)
          ijcov(i,1)=j
          ijcov(i,2)=j
          j=j+1
        enddo

c
c  Add smoothing equations
        if (rank.eq.0) write(*,*) 'Define regularisation'
        call build_damp_space(nlocdatpts, nlocsampts,
     >                        imin_locpts, imin_locsampts,
     >                        ND, ncoeffs, shdeg, dampfac,
     >                        bc, ijcov, cov, ppos)

c
c  Prepare the CM4 components for use within sub_sph_wmam_l()
        do i=1,nlocpts
          call prepare_cm4_components(ppos(1,i))
        enddo
        
c
c  Finalise covariance matrix
        call DS2Y(nlocpts,nlocpts,ijcov(1,1),ijcov(1,2),cov,0)

c
c  Read in starting model
        fname='./Data/model.in.bin'
        if (rank.eq.0) write(*,*)
     >    'Reading in starting model, ', fname
        if (serialrd .gt. 1) then
          call mpi_read_ini_model(fname, nparams, bc)
        else
          call mpi_read_all_ini_model(fname, nparams, bc)
        endif
        if (nparams .eq. 0) stop

c
c  Invert data
        itmax(1)=7
        itmax(2)=5
        itmax(3)=-10
        stdt=1.0d0
        dl(1)=1.0d-10
        dl(2)=0.0d0
        dl(3)=1.d14
c
        dw=1.d0
c 
        if (rank.eq.0) then
          write(*,*) 'Start Inversion'
          write(*,*) ' itmax: ', itmax(1:3)
          write(*,*) ' dl: ', dl(1:3)
          write(*,*) ''
         endif
c
        fname='./Results/'
        allocate(gg(1:1,1:1))
        allocate(bb(1:1))

c
        if (scheme.eq.POLAK_RIBIERE) then
          call opt_pr_p3(fname, itmax, NPMAX, ND, nparams,
     >                   npts, nlocpts, ppos, bc, dl,
     >                   sub_base_i, damp_rien,
     >                   fun_base_f, l2_norm, l2_std,
     >                   cov, ijcov(1,2),
     >                   stdt, dw, bb, gg)
        else
c         CONJUGATE_GRADIENT
          call opt_ghc_p2(fname, itmax, NPMAX, ND, nparams,
     >                    npts, nlocpts, ppos, bc, dl,
     >                    sub_base_i, damp_rien,
     >                    fun_base_f, l2_norm, l2_std,
     >                    cov, ijcov(1,2),
     >                    stdt, dw, bb, gg)
        endif

c
        deallocate(bb)
        deallocate(gg)
c
        if (rank.eq.0) then
          write(*,'(A)')' '
          write(*,'(A,e15.7)') 'The L2 STD is: ',stdt
          write(*,'(A)')' '
        endif
c
c
        fname='./Results/fit_No_P.out.bin'
        call mpi_write_all_fit_data(fname, ND,
     >                              ndatpts, imin_locdatpts,
     >                              1, nlocdatpts,
     >                              ppos, dw, diff(1:2))

        fname='./Results/fit_damp.out.bin'
        call mpi_write_all_fit_data(fname, ND,
     >                              nsampts, imin_locsampts,
     >                              nlocdatpts+1, nlocpts,
     >                              ppos, dw, diff(3:4))
c
        if (rank .eq. 0) then
          call MPI_Reduce(MPI_IN_PLACE, diff, 4,
     >                    MPI_DOUBLE_PRECISION,
     >                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        else
          call MPI_Reduce(diff, diff, 4,
     >                    MPI_DOUBLE_PRECISION,
     >                    MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        endif
c
c       
        if (rank .eq. 0) then
          diff(1)=diff(1)/ndatpts
          diff(2)=dsqrt((diff(2)-ndatpts*diff(1)**2)/ndatpts)
c
          write(*,'(2(A,f16.7))')'Residual average : ',diff(1)
     >                          ,' with L2 std : ',diff(2)

          diff(3)=diff(3)/nsampts
          diff(4)=dsqrt((diff(4)-nsampts*diff(3)**2)/nsampts)
c
          write(*,'(2(A,f16.7))')'Residual average : ',diff(3)
     >                          ,' with L2 std : ',diff(4)

c
          std=stdt*npts/(npts-nparams)
          allocate(err(1:nparams))
          err(1:nparams)=std
c
c  Saving update base coefficients
          fname='./Results/model_No_P.out'
          open(10,file=fname)
            write(10,'(A)') '#'
            write(10,'(A,i8)') '#lmax= ', shdeg
            write(10,'(A)') '###########'
            write(10,*)ryg
            i=0
            do il=1,shdeg
              im=0
              i=i+1
              write(10,*) 'c', il, im, bc(i), err(i)
              do im=1,il
                i=i+1
                write(10,*) 'c', il, im, bc(i), err(i)
                i=i+1
                write(10,*) 'c', il, -im, bc(i), err(i)
              enddo
            enddo
          close(10)

          deallocate(err)
        endif
c
c  Deallocate arrays
        deallocate(dw)
        deallocate(ijcov)
        deallocate(cov)
        deallocate(ppos)
        deallocate(bc)
        deallocate(proc_np,proc_ip)
c
        call fini_sph_wmam()
c
        call MPI_Finalize(ierr)
c
        end
