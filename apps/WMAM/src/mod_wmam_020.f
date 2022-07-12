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
        character(len=*), parameter :: VERSION="2.2.0"
c
        integer, parameter :: POLAK_RIBIERE=1
        integer, parameter :: CONJUGATE_GRADIENT=2
        integer, parameter :: ND=7
        integer, parameter :: NPMAX=3
c
        character fname*100
        character buf*100
        integer i, j, ix, iy, il, im
        integer nx, ny
        integer itmax(3), nub
        real*8  std, stdt, dd, dl(3), wgh
        real*8  resdeg
        real*8  longitude, colatitude
        real*8  radius, magscalar
        integer shdeg, scheme
        integer ncoeffs, nparams
        integer ndatpts, nsampts, npts
        integer nlocdatpts, imin_locdatpts, imax_locdatpts
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
        if (COMMAND_ARGUMENT_COUNT().eq.3) then
          call GET_COMMAND_ARGUMENT(1,argstr)
          read(argstr,*) shdeg
          call GET_COMMAND_ARGUMENT(2,argstr)
          read(argstr,*) resdeg
          resdeg = nint(resdeg * 1000.0) * 1E-3
          call GET_COMMAND_ARGUMENT(3,argstr)
          read(argstr,*) scheme
        else
          shdeg=200
          resdeg=1.0
          scheme=POLAK_RIBIERE
        endif
        if (rank.eq.0) then
          write(*,*) 'WMAM v', VERSION
          write(*,*) ''
          write(*,*) 'shdeg: ', shdeg
          write(*,*) 'resdeg: ', resdeg
          write(*,*) 'scheme: ', scheme
          write(*,*) ''
        endif
c
c  Settings
        call init_sph_wmam(shdeg)
        nparams=shdeg*(shdeg+2)
        ny=int(1.0/resdeg)*360
        nx=int(1.0/resdeg)*180-1
        ndatpts=nx*ny
        nsampts=(shdeg+1)*(2*shdeg+1)
        ncoeffs=nparams
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
        imax_locdatpts = imin_locdatpts + nlocdatpts - 1

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
c  Array allocations
        allocate(bc(nparams))
        allocate(ppos(ND+1,nlocpts))
        allocate(cov(nlocpts))
        allocate(ijcov(nlocpts,2))
        allocate(dw(nlocpts))
        
c
c  Read in reference field
        bc(1:nparams)=1.0d0
        fname='./Data/coef_1990_15.dat'
        if (rank.eq.0) write(*,*)
     >    'Reading in reference field, ', fname
        call read_model(fname, ryg, bc, ncoeffs)
        if (rank.eq.0) then
          write(*,*) 'Coefficients: ', ncoeffs
          write(*,*) ''
        endif

c  Read in data
        fname='./Data/wdmam_geocentric.dat'
        if (rank.eq.0) then
          write(*,*) 'Reading in data, ', fname
          write(*,*) ''
        endif
        open(10,file=fname,status='old')
          i=1
          j=1
          do ix=1,nx
            do iy=1,ny
              read(10,*) longitude, colatitude, radius, magscalar

              if (i .ge. imin_locdatpts .and.
     >            i .le. imax_locdatpts) then
                ppos(1,j) = 90.0d0 - colatitude
                ppos(2,j) = longitude
                ppos(3,j) = radius
                ppos(4,j) = ryg
                ppos(ND+1,j) = magscalar
                j=j+1
              endif

              i=i+1
            enddo
          enddo
        close(10)
        
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
        wgh=5.0d0
        call build_damp_space(nlocdatpts, nlocsampts,
     >                        imin_locpts, imin_locsampts,
     >                        ND, ncoeffs, shdeg, wgh, bc,
     >                        ijcov, cov, ppos)
c

c
c  Finalise covariance matrix
        call DS2Y(nlocpts,nlocpts,ijcov(1,1),ijcov(1,2),cov,0)

c
c  Read in starting model
        fname='./Data/model.in'
        if (rank.eq.0) then
          write(*,*) 'Reading in starting model, ', fname
          write(*,*) ''
        endif
        bc(1:nparams)=0.0d0
        open(10,file=fname,status='old')
          read(10,*)buf
          read(10,*)buf
          read(10,*)buf
          read(10,*)buf
          do i=1,nparams
            read(10,*) j, bc(i), dd
          enddo
        close(10)

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
          write(*,*) ' npts: ', npts
          write(*,*) ' wgh: ', wgh
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
     >                   nlocdatpts, proc_np, ppos, bc, dl,
     >                   l2_norm, sub_sph_wmam_l, l2_std, damp_rien,
     >                   cov, ijcov(1,2),
     >                   stdt, dw, bb, gg)
        else
c         CONJUGATE_GRADIENT
          call opt_ghc_p2(fname, itmax, NPMAX, ND, nparams,
     >                    nlocdatpts, proc_np, ppos, bc, dl,
     >                    l2_norm, sub_sph_wmam_l, l2_std, damp_rien,
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
c
c  Writing fit to data per component
          fname='./Results/fit_No_P.out'
          call write_comp(fname,nlocpts,ND,1,nlocdatpts,ppos,dw)
          fname='./Results/fit_damp.out'
          call write_comp(fname,nlocpts,ND,nlocdatpts+1,nlocpts,
     >                    ppos,dw)
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
        deallocate(ppos)
        deallocate(bc)
        deallocate(cov)
        deallocate(ijcov)
        deallocate(dw)
        deallocate(proc_np,proc_ip)
c
        call fini_sph_wmam()
c
        call MPI_Finalize(ierr)
c
        stop
        end
