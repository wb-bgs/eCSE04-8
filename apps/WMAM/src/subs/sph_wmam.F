cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     module sph_wmam
c                       V. Lesur 16/09/2006
c
c     This module set the variable and define the model for the
c     anaomaly field data
c
c     All data are scalar
c
c  Number of parameters as organised on BE:
c  ilg*(ilg+2)                        internal lithosphere models
c
c    Global data:
c         ilg : integer : internal SH degree value (see below)
c  imax_datpt : integer : number of local data points
c          nc : integer : number of points in cache
c         ryg : real*8  : reference year for the model
c         rag : real*8  : radius for lithosphere
c         dwx, dwy, dwz : real*8 : allocatable arrays of size nb,0:ncache
c         d2a : real*8  : allocatable array of size 0:ilg;
c                         holds pre-computed values for mklf_F2()
c    Subroutines:
c       init_sph_wmam()    initialise global variables
c       fini_sph_wmam()    deallocate global arrays
c       mpi_read_ref_model()
c       mpi_read_model()
c       mpi_read_data()
c       prepare_cm4_components()
c       populate_dw_cache()
c       sub_sph_wmam_l()   linearized model subroutine
c       sub_base_i()
c    Functions
c       fun_base_f()
c       wmam_norm()
c       wmam_var
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        module sph_wmam
c
        implicit none
c
        integer, parameter :: NUM_OF_PTCOMPS = 4
        integer, parameter :: SIZE_OF_REAL = 8
        integer, parameter :: SIZE_OF_PT = NUM_OF_PTCOMPS*SIZE_OF_REAL
c        
        integer ilg, imax_datpt, nc
        real*8 ryg, rag
c
        real*8, allocatable :: dwx(:,:), dwy(:,:), dwz(:,:)
        real*8, allocatable :: d2a(:)
c
        contains
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine init_sph_wmam
c
c       V. Lesur 16/09/2006
c
c  Set ilg, ryg (based on 1990) and rag
c
c  input:
c       shdeg : integer : spherical harmonic degree
c          nb : integer : number of parameters
c         nlp : integer : number of local data points
c      ncache : integer : number of points in cache
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine init_sph_wmam(shdeg, nb, nlp, ncache)
c
        implicit none
c
        integer shdeg, nb, nlp, ncache
        integer nm, ierr
        real*8  dnm, d1, d2
        real*8 dgamln
c
        external dgamln 
c
        call dy2mjd(1990.d0,ryg)
        rag=6371.2d0
        ilg=shdeg
c
        allocate(d2a(0:ilg))
c
        do nm=0,ilg
          dnm = dble(nm)                    ! dble real for nm
          d1 = dgamln(2*dnm+1.0d0,ierr)     ! d1=log(fact(2dnm))
          d2 = dgamln(dnm+1.0d0,ierr)       ! d2=log(fact(dnm))
          if (ierr.ne.0) then
            write(*,*) 'mklf_F: Cannot computes normalisation cst !'
            stop
          endif
c
          d2 = 0.5d0*d1 - d2                ! d2=sqrt(fact(2dnm))/fact(dnm)
          d2 = d2 - nm*dlog(2.0d0)          !
          d2 = dexp(d2)                     ! normalisation cst.
          if (nm.ne.0) d2 = d2*dsqrt(2.0d0) ! special case  m=0
c
          d2a(nm) = d2
        enddo
c
        imax_datpt = nlp
        nc = ncache
c  allocate cache 
        allocate(dwx(nb, 0:nc))
        allocate(dwy(nb, 0:nc))
        allocate(dwz(nb, 0:nc))
c
        return
        end subroutine init_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine fini_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fini_sph_wmam()
        deallocate(dwx,dwy,dwz)
        deallocate(d2a)
c
        nc = 0
c
        return
        end subroutine fini_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_ref_model
c
c  Read in reference model.
c  All ranks read entire file.
c
c  input:
c       fname   : character : name of reference model file
c  output:
c       ncoeffs : integer : number of coefficients/parameters
c       bc(*)   : real*8  : coefficients/parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_ref_model(fname, ncoeffs, bc)
c
        implicit none
c
        include 'mpif.h'
c
	      character fname*100
        integer ncoeffs
        real*8 bc(*)
c
        integer, parameter :: ndims = 1
c
        integer, dimension(ndims) :: array_of_sizes
        integer, dimension(ndims) :: array_of_subsizes
        integer, dimension(ndims) :: array_of_starts
        integer subarray
c
        integer rank, ierr
        integer fhandle, fsize, nreals, nread
        integer (kind=MPI_OFFSET_KIND) :: disp = 0
        integer fstat(MPI_STATUS_SIZE)
c
        logical file_open, buf_allocated
        logical subarray_type_committed
c       
        character mpifunc*40 
        real*8, allocatable :: buf(:)
c
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
        file_open = .false.
        buf_allocated = .false.
        subarray_type_committed = .false.
c
        do while (.true.)
c
          call MPI_File_open(MPI_COMM_WORLD, fname,
     >                       MPI_MODE_RDONLY, MPI_INFO_NULL,
     >                       fhandle, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            file_open = .true.
          else
            mpifunc = 'MPI_File_open'
            exit
          endif 
c        
          call MPI_File_get_size(fhandle, fsize, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            nreals = fsize / SIZE_OF_REAL
          else
            nreals = 0
            mpifunc = 'MPI_File_get_size'
            exit
          endif

          if (nreals .gt. 1+ncoeffs) then
            if (rank .eq. 0) then
              write(*,*) 'Error, ref model file should have ',
     >                   'no more than 1 +', ncoeffs, ' reals.'
            endif
            ncoeffs = 0
            exit
          endif

          allocate(buf(nreals))
          buf_allocated = .true.
          
          array_of_sizes(1) = nreals
          array_of_subsizes(1) = array_of_sizes(1)
          array_of_starts(1) = 0

          call MPI_Type_create_subarray(ndims, array_of_sizes,
     >                                  array_of_subsizes,
     >                                  array_of_starts,
     >                                  MPI_ORDER_FORTRAN,
     >                                  MPI_DOUBLE_PRECISION,
     >                                  subarray, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_Type_create_subarray'
            exit
          endif

          call MPI_Type_commit(subarray, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            subarray_type_committed = .true.
          else
            mpifunc = 'MPI_Type_commit'
            exit 
          endif
          
          call MPI_File_set_view(fhandle, disp,
     >                           MPI_DOUBLE_PRECISION,
     >                           subarray, 'native', 
     >                           MPI_INFO_NULL, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_set_view'
            exit
          endif

          call MPI_File_read_all(fhandle, buf, 1,
     >                           subarray, fstat, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_read_all'
            exit
          endif

          call MPI_Get_count(fstat, subarray, nread, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            if (nread .eq. 1) then
              call dy2mjd(buf(1), ryg)
              ncoeffs = nreals-1
              bc(1:ncoeffs) = buf(2:nreals)
            else
              if (rank .eq. 0) then
                write(*,*) 'Error, unable to read 1 +', ncoeffs,
     >                     ' reals from ref model file.'
              endif
              ncoeffs = 0
            endif
          else
            mpifunc = 'MPI_Get_count'
            exit
          endif

          exit
        enddo

        if (ierr .ne. MPI_SUCCESS) then
          if (rank .eq. 0) then
            write(*,*) 'Error ', mpifunc, '() returned ', ierr
          endif
          ncoeffs = 0
        endif

        if (subarray_type_committed) then
          call MPI_Type_free(subarray, ierr)
        endif

        if (buf_allocated) then
          deallocate(buf)
        endif

        if (file_open) then
          call MPI_File_close(fhandle, ierr)
        endif

        end subroutine mpi_read_ref_model
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_model
c
c  Read in starting model.
c  All ranks read entire file.
c
c  input:
c       fname   : character : name of model file
c  output:
c       nparams : integer : number of parameters
c       bc(*)   : real*8  : parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_model(fname, nparams, bc)
c
        implicit none
c
        include 'mpif.h'
c
	      character fname*100
        integer nparams
        real*8 bc(*)
c
        integer, parameter :: ndims = 1
c
        integer, dimension(ndims) :: array_of_sizes
        integer, dimension(ndims) :: array_of_subsizes
        integer, dimension(ndims) :: array_of_starts
        integer subarray
c
        integer rank, ierr
        integer fhandle, fsize, nreals, nread
        integer (kind=MPI_OFFSET_KIND) :: disp = 0
        integer fstat(MPI_STATUS_SIZE)
c
        logical file_open
        logical subarray_type_committed
c       
        character mpifunc*40  
c
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
        file_open = .false.
        subarray_type_committed = .false.
c
        do while (.true.)
c
          call MPI_File_open(MPI_COMM_WORLD, fname,
     >                       MPI_MODE_RDONLY, MPI_INFO_NULL,
     >                       fhandle, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            file_open = .true.
          else
            mpifunc = 'MPI_File_open'
            exit
          endif 
c        
          call MPI_File_get_size(fhandle, fsize, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            nreals = fsize / SIZE_OF_REAL
          else
            nreals = 0
            mpifunc = 'MPI_File_get_size'
            exit
          endif

          if (nreals .ne. nparams) then
            if (rank .eq. 0) then
              write(*,*) 'Error, model file should have',
     >                   nparams, ' reals.'
            endif
            nparams = 0
            exit
          endif

          array_of_sizes(1) = nparams
          array_of_subsizes(1) = array_of_sizes(1)
          array_of_starts(1) = 0

          call MPI_Type_create_subarray(ndims, array_of_sizes,
     >                                  array_of_subsizes,
     >                                  array_of_starts,
     >                                  MPI_ORDER_FORTRAN,
     >                                  MPI_DOUBLE_PRECISION,
     >                                  subarray, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_Type_create_subarray'
            exit
          endif

          call MPI_Type_commit(subarray, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            subarray_type_committed = .true.
          else
            mpifunc = 'MPI_Type_commit'
            exit 
          endif
          
          call MPI_File_set_view(fhandle, disp,
     >                           MPI_DOUBLE_PRECISION,
     >                           subarray, 'native', 
     >                           MPI_INFO_NULL, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_set_view'
            exit
          endif

          call MPI_File_read_all(fhandle, bc, 1,
     >                           subarray, fstat, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_read_all'
            exit
          endif

          call MPI_Get_count(fstat, subarray, nread, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            if (nread .ne. 1) then
              if (rank .eq. 0) then
                write(*,*) 'Error, unable to read', nparams,
     >                     ' reals from ref model file.'
              endif
              nparams = 0
            endif
          else
            mpifunc = 'MPI_Get_count'
            exit
          endif

          exit
        enddo

        if (ierr .ne. MPI_SUCCESS) then
          if (rank .eq. 0) then
            write(*,*) 'Error ', mpifunc, '() returned ', ierr
          endif
          nparams = 0
        endif

        if (subarray_type_committed) then
          call MPI_Type_free(subarray, ierr)
        endif

        if (file_open) then
          call MPI_File_close(fhandle, ierr)
        endif

        end subroutine mpi_read_model
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_data
c
c  Read in data.
c  Read just those points that have been assigned to the rank.
c
c  input:
c       fname      : character : name of data file
c       nd         : integer : nd+1 is the lead dim of ppos
c       ndatpts    : integer : number of data points in file
c       nlocdatpts : integer : number of data points assigned to rank
c       imin_locdatpts : integer : one-based rank index for global
c                                  array of data points
c  output:
c       ppos(*,*) : real*8 : array of data points assigned to rank
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_data(fname, nd, ndatpts, nlocdatpts,
     >                           imin_locdatpts, ppos)
c
        implicit none
c
        include 'mpif.h'
c
	      character fname*100
        integer nd, ndatpts, nlocdatpts, imin_locdatpts
        real*8 ppos(nd+1,*)
c
        integer i, j, rank, ierr
        integer fhandle, fsize, npts, nread, nlocreals
        integer (kind=MPI_OFFSET_KIND) :: disp = 0
        integer fstat(MPI_STATUS_SIZE)
c
        logical file_open, buf_allocated
c       
        character mpifunc*40 
        real*8, allocatable :: buf(:)
c
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
        file_open = .false.
        buf_allocated = .false.
c
        do while (.true.)
c
          call MPI_File_open(MPI_COMM_WORLD, fname,
     >                       MPI_MODE_RDONLY, MPI_INFO_NULL,
     >                       fhandle, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            file_open = .true.
          else
            mpifunc = 'MPI_File_open'
            exit
          endif
c        
          call MPI_File_get_size(fhandle, fsize, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            npts = (fsize / SIZE_OF_PT)
          else
            npts = 0
            mpifunc = 'MPI_File_get_size'
            exit
          endif

          if (npts .ne. ndatpts) then
            if (rank .eq. 0) then
              write(*,*) 'Error, data file should have',
     >                   ndatpts, ' points.'
            endif
            nlocdatpts = 0
            exit
          endif

          nlocreals = nlocdatpts*NUM_OF_PTCOMPS
          allocate(buf(nlocreals))
          buf_allocated = .true.
 
          disp = (imin_locdatpts-1)*SIZE_OF_PT

          call MPI_File_set_view(fhandle, disp,
     >                           MPI_DOUBLE_PRECISION,
     >                           MPI_DOUBLE_PRECISION,
     >                           'native', 
     >                           MPI_INFO_NULL, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_set_view'
            exit 
          endif

          call MPI_File_read_all(fhandle, buf, nlocreals,
     >                           MPI_DOUBLE_PRECISION,
     >                           fstat, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_read_all'
            exit
          endif

          call MPI_Get_count(fstat, MPI_DOUBLE_PRECISION, 
     >                       nread, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            if (nread .eq. nlocreals) then
              j=1
              do i=1,nlocdatpts
                ppos(1,i) = 90.0d0 - buf(j+1)
                ppos(2,i) = buf(j)
                ppos(3,i) = buf(j+2)
                ppos(4,i) = ryg
                ppos(nd+1,i) = buf(j+3)
                j=j+4
              enddo
            else
              if (rank .eq. 0) then
                write(*,*) 'Error, unable to read', nlocdatpts,
     >                     ' points from data file.'
              endif
              nlocdatpts = 0
            endif
          else
            mpifunc = 'MPI_Get_count'
            exit
          endif

          exit
        enddo

        if (ierr .ne. MPI_SUCCESS) then
          if (rank .eq. 0) then
            write(*,*) 'Error ', mpifunc, '() returned ', ierr
          endif
          nlocdatpts = 0
        endif

        if (buf_allocated) then
          deallocate(buf)
        endif

        if (file_open) then
          call MPI_File_close(fhandle, ierr)
        endif

        end subroutine mpi_read_data
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine prepare_cm4_components
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine prepare_cm4_components(bp)

       implicit none

       real*8 bp(*)
       real*8 dd

       dd = dsqrt(bp(5)**2 + bp(6)**2 + bp(7)**2)
       bp(5) = bp(5) / dd
       bp(6) = bp(6) / dd
       bp(7) = bp(7) / dd

       end subroutine prepare_cm4_components
c

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine populate_dw_cache
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine populate_dw_cache(nd, ppos)

       implicit none

       integer nd
       real*8  ppos(nd+1,*)

       integer i

       do i=1,nc
         call XYZsph_bi0(ilg, rag, ppos(1,i), d2a,
     >                   dwx(1,i), dwy(1,i), dwz(1,i))
       enddo

       end subroutine populate_dw_cache
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sub_sph_wmam_l
c
c       V. Lesur  16/09/2006
c
c   That is for a linearized inversion of the lithosphere field.
c   This subroutine is called from either sub_base_i() or fun_base_f().
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine sub_sph_wmam_l(ipt, nb, bc, bp, be)
c
       implicit none
c
#ifdef CRAYPAT
       include 'pat_apif.h'
#endif
c
       integer ipt, nb
       real*8 bc(*), bp(*), be(*)
c
       integer i
c
       real*8 dx, dy, dz, dd
c
       real*8 dxbey, dxbez
       real*8 dybex, dybez
       real*8 dzbex, dzbey
c
       real*8 xy_c, xz_c
       real*8 yx_c, yz_c
       real*8 zx_c, zy_c
c
       real*8 bex, bey, bez
       real*8 bex2, bey2, bez2

#ifdef CRAYPAT
       integer cpstat

       call PAT_record(PAT_STATE_ON,cpstat)
       call PAT_region_begin(1,"XYZsph_bi0",cpstat)
#endif

c    calculate internal field component
       if (ipt .gt. nc) then
         call XYZsph_bi0(ilg, rag, bp, d2a,
     >                   dwx(1,0), dwy(1,0), dwz(1,0))
         i = 0
       else
c    dw arrays have been pre-calculated for this point
         i = ipt
       endif

#ifdef CRAYPAT
       call PAT_region_end(1,cpstat)
       call PAT_record(PAT_STATE_OFF,cpstat)
#endif

       bex = bp(5)
       bey = bp(6)
       bez = bp(7)

       if (ipt .gt. imax_datpt) then
c    sampling point
         dx = dot_product(dwx(1:nb,i), bc(1:nb))
         dy = dot_product(dwy(1:nb,i), bc(1:nb))
         dz = dot_product(dwz(1:nb,i), bc(1:nb))
c
         dxbey = dx*bey
         dxbez = dx*bez
         dybex = dy*bex
         dybez = dy*bez
         dzbex = dz*bex
         dzbey = dz*bey
c
         xy_c = dxbey - dybex
         xz_c = dxbez - dzbex
         yx_c = -xy_c
         yz_c = dybez - dzbey
         zx_c = -xz_c
         zy_c = -yz_c
c
         dd = dsqrt(yz_c**2 + xz_c**2 + xy_c**2)
c
         bex2 = (xz_c*bez + xy_c*bey) / dd
         bey2 = (yz_c*bez + yx_c*bex) / dd
         bez2 = (zy_c*bey + zx_c*bex) / dd

         bex = bex2
         bey = bey2
         bez = bez2
       endif
c
       be(1:nb) = dwx(1:nb,i)*bex
     >          + dwy(1:nb,i)*bey
     >          + dwz(1:nb,i)*bez
c
       return
       end subroutine sub_sph_wmam_l
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sub_base_i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine sub_base_i(ipt, nb, bc, bp, be)
c
       implicit none
c
       integer ipt, nb
       real*8 bc(*), bp(*), be(*)
c
       call sub_sph_wmam_l(ipt, nb, bc, bp, be)
c 
       return
       end subroutine sub_base_i
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    function fun_base_f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 function fun_base_f(ipt, nb, bc, bp, be)
c
       implicit none
c
       integer ipt, nb
       real*8 bc(*), bp(*), be(*)
c
       call sub_sph_wmam_l(ipt, nb, bc, bp, be)
c
       fun_base_f = dot_product(be(1:nb), bc(1:nb))
c
       return
       end function fun_base_f
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 function wmam_norm(i,nub,mv)
c
       implicit none
c
       integer i,nub(*)
       real*8 mv(*)
c
       real*8 dc,ae
c
       real*8 lesur_norm,l2_norm,l1_norm
       external lesur_norm,l2_norm,l1_norm
c
       select case (nub(i))
       case (1)
         dc=0.6d0
         ae=0.5d0
         wmam_norm=l2_norm(i,nub,mv)
c        wmam_norm=lesur_norm(dc,ae,i,mv)
       case default
         wmam_norm=l2_norm(i,nub,mv)
c        wmam_norm=l1_norm(i,nub,mv)
       end select
c
       return
       end function wmam_norm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       real*8 function wmam_var(nub,npt,np,var,mv,wgh)
c
       implicit none
c
       integer npt,np,nub(*)
       real*8 mv(*),var,wgh(*)
c
       integer i,j,k
       real*8 vv,dv
       real*8 dc,ae
c
       real*8 lesur_var2,l2_var,l1_var
       external lesur_var2,l2_var,l1_var
c
       if(np.ge.0)then
         vv=var
         do i=1,np
           dv=vv
           j=0
           k=1
c
           select case (nub(i))
           case (1)
             dc=0.6d0
             ae=0.5d0
             vv=l2_var(nub(i),j,k,dv,mv(i),wgh(i))
c            vv=lesur_var2(dc,ae,j,k,dv,mv(i),wgh(i))
           case default
             vv=l2_var(nub(i),j,k,dv,mv(i),wgh(i))
c            vv=l1_var(nub(i),j,k,dv,mv(i),wgh(i))
           end select
         enddo
       else
         vv=0.0d0
         do i=1,-np
           vv=vv+mv(i)
         enddo
       endif
c
       wmam_var=vv
c
       return
       end function wmam_var
ccccccccccccccccccccccccccccccccccc
       end module sph_wmam