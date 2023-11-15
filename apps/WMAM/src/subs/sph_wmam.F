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
c        ilg : integer : internal SH degree value (see below)
c       nldp : integer : number of data points assigned to rank
c        ryg : real*8  : reference year for the model
c        rag : real*8  : radius for lithosphere
c       d2a : real*8  : allocatable array of size 0:ilg;
c                       holds pre-computed values for mklf_F2()
c    Subroutines:
c       init_sph_wmam()    initialise global variables
c       fini_sph_wmam()    deallocate global arrays
c       mpi_read_ref_model()
c       mpi_read_ini_model()
c       mpi_read_all_ref_model()
c       mpi_read_all_ini_model()
c       mpi_read_all_data()
c       mpi_write_all_fit_data()
c       prepare_cm4_components()
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
        integer, parameter :: FIT_VALUES_PER_PT = 10
c        
        integer ilg, nldp
        real*8 ryg, rag
c
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
c             shdeg : integer : spherical harmonic degree
c                nb : integer : number of parameters
c        nlocdatpts : integer : number of data points assigned to rank
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine init_sph_wmam(shdeg, nb, nlocdatpts)
c
        implicit none
c
        integer shdeg, nb, nlocdatpts
        integer nm, ierr
        real*8  dnm, d1, d2
        real*8 dgamln
c
        external dgamln 
c
        rag=6371.2d0
        ilg=shdeg
        nldp=nlocdatpts
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
        return
        end subroutine init_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine fini_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fini_sph_wmam()
        deallocate(d2a)
        return
        end subroutine fini_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_ref_model
c
c  Read in reference model.
c  Rank 0 does the read and then broadcasts to all other ranks.
c
c  input:
c       fname   : character : name of reference model file
c       ncoeffs : integer : number of coefficients/parameters
c  output:
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
        integer rank, ierr
        integer fhandle, fsize, nreals, nread
        integer (kind=MPI_OFFSET_KIND) :: disp = 0
        integer fstat(MPI_STATUS_SIZE)
c
        logical file_open
c       
        character mpifunc*40 
        real*8, allocatable :: buf(:)
c
c
        nreals = 1 + ncoeffs
        allocate(buf(nreals))
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
        if (rank .eq. 0) then
c
          file_open = .false.
c
          do while (.true.)
c
            call MPI_File_open(MPI_COMM_SELF, fname,
     >                         MPI_MODE_RDONLY, MPI_INFO_NULL,
     >                         fhandle, ierr)
            if (ierr .eq. MPI_SUCCESS) then
              file_open = .true.
            else
              mpifunc = 'MPI_File_open'
              exit
            endif 
c        
            call MPI_File_get_size(fhandle, fsize, ierr)
            if (ierr .eq. MPI_SUCCESS) then
              nread = fsize / SIZE_OF_REAL
            else
              mpifunc = 'MPI_File_get_size'
              exit
            endif
c
            if (nread .ne. nreals) then
              write(*,*) 'Error, ref model file should contain ',
     >                   '1 +', ncoeffs, ' reals.'
              nread = 0
              exit
            endif
c
            call MPI_File_set_view(fhandle, disp,
     >                             MPI_DOUBLE_PRECISION,
     >                             MPI_DOUBLE_PRECISION,
     >                             'native', 
     >                             MPI_INFO_NULL, ierr)
            if (ierr .ne. MPI_SUCCESS) then
              mpifunc = 'MPI_File_set_view'
              exit
            endif
c
            call MPI_File_read(fhandle, buf, nreals,
     >                         MPI_DOUBLE_PRECISION,
     >                         fstat, ierr)
            if (ierr .ne. MPI_SUCCESS) then
              mpifunc = 'MPI_File_read_all'
              exit
            endif
c
            call MPI_Get_count(fstat, MPI_DOUBLE_PRECISION,
     >                         nread, ierr)
            if (ierr .eq. MPI_SUCCESS) then
              if (nread .ne. nreals) then
                write(*,*) 'Error, unable to read 1 +', ncoeffs,
     >                     ' reals from ref model file.'
              endif
            else
              mpifunc = 'MPI_Get_count'
              exit
            endif
c
            exit
c
          enddo
c
          if (file_open) then
            call MPI_File_close(fhandle, ierr)
          endif
c
          if (ierr .eq. MPI_SUCCESS) then
            if (nread .eq. nreals) then
              call MPI_BCAST(buf, nreals, MPI_DOUBLE_PRECISION,
     >                       0, MPI_COMM_WORLD, ierr)
            else
              stop
            endif
          else
            write(*,*) 'Error ', mpifunc, '() returned ', ierr
            stop
          endif
c          
        else
c  Rank is not zero
          call MPI_BCAST(buf, nreals, MPI_DOUBLE_PRECISION,
     >                   0, MPI_COMM_WORLD, ierr)
        endif
c
        call dy2mjd(buf(1), ryg)
        bc(1:ncoeffs) = buf(2:nreals)
c
        deallocate(buf)
c
        end subroutine mpi_read_ref_model
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_ini_model
c
c  Read in starting model.
c  Rank zero reads in file and then broadcasts to all other ranks.
c
c  input:
c       fname   : character : name of model file
c       nparams : integer : number of parameters
c  output:
c       bc(*)   : real*8  : parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_ini_model(fname, nparams, bc)
c
        implicit none
c
        include 'mpif.h'
c
        character fname*100
        integer nparams
        real*8 bc(*)
c
        integer rank, ierr
        integer fhandle, fsize, nreals, nread
        integer (kind=MPI_OFFSET_KIND) :: disp = 0
        integer fstat(MPI_STATUS_SIZE)
c
        logical file_open
c       
        character mpifunc*40  
c
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
        if (rank .eq. 0) then
c
          file_open = .false.
c
          do while (.true.)
c
            call MPI_File_open(MPI_COMM_SELF, fname,
     >                         MPI_MODE_RDONLY, MPI_INFO_NULL,
     >                         fhandle, ierr)
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
              mpifunc = 'MPI_File_get_size'
              exit
            endif
c
            if (nreals .ne. nparams) then
              write(*,*) 'Error, model file should have',
     >                   nparams, ' reals.'
              nread = 0
              exit
            endif
c
            call MPI_File_set_view(fhandle, disp,
     >                             MPI_DOUBLE_PRECISION,
     >                             MPI_DOUBLE_PRECISION,
     >                             'native', 
     >                             MPI_INFO_NULL, ierr)
            if (ierr .ne. MPI_SUCCESS) then
              mpifunc = 'MPI_File_set_view'
              exit
            endif
c
            call MPI_File_read(fhandle, bc, nparams,
     >                         MPI_DOUBLE_PRECISION,
     >                         fstat, ierr)
            if (ierr .ne. MPI_SUCCESS) then
              mpifunc = 'MPI_File_read_all'
              exit
            endif
c
            call MPI_Get_count(fstat, MPI_DOUBLE_PRECISION,
     >                         nread, ierr)
            if (ierr .eq. MPI_SUCCESS) then
              if (nread .ne. nparams) then
                write(*,*) 'Error, unable to read', nparams,
     >                     ' reals from ref model file.'
              endif
            else
              mpifunc = 'MPI_Get_count'
              exit
            endif
c
            exit
c
          enddo
c
          if (file_open) then
            call MPI_File_close(fhandle, ierr)
          endif
c
          if (ierr .eq. MPI_SUCCESS) then
            if (nread .eq. nparams) then
              call MPI_BCAST(bc, nparams, MPI_DOUBLE_PRECISION,
     >                       0, MPI_COMM_WORLD, ierr)
            else
              stop
            endif
          else
            write(*,*) 'Error ', mpifunc, '() returned ', ierr
            stop
          endif
c         
        else
c  Rank is not zero
          call MPI_BCAST(bc, nparams, MPI_DOUBLE_PRECISION,
     >                   0, MPI_COMM_WORLD, ierr)
        endif
c
        end subroutine mpi_read_ini_model
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_all_ref_model
c
c  Read in reference model.
c  All ranks read entire file via MPI_File_read_all().
c
c  input:
c       fname   : character : name of reference model file
c  output:
c       ncoeffs : integer : number of coefficients/parameters
c       bc(*)   : real*8  : coefficients/parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_all_ref_model(fname, ncoeffs, bc)
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

        end subroutine mpi_read_all_ref_model
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_all_ini_model
c
c  Read in starting model.
c  All ranks read entire file via MPI_File_read_all().
c
c  input:
c       fname   : character : name of model file
c  output:
c       nparams : integer : number of parameters
c       bc(*)   : real*8  : parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_all_ini_model(fname, nparams, bc)
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

        end subroutine mpi_read_all_ini_model
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_all_data
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
        subroutine mpi_read_all_data(fname, nd, ndatpts, nlocdatpts,
     >                               imin_locdatpts, ppos)
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

        end subroutine mpi_read_all_data
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_write_all_fit_data
c
c  Write out fit data.
c  Write out the values for just those points that have been assigned to the rank.
c
c  input:
c       fname      : character : name of fit data file
c       nd         : integer : nd+1 is the lead dim of ppos
c       np         : integer : total number of points
c       imin       : integer : offset into fit data file.
c       nmin       : integer : starting count
c       nmax       : integer : finishing count
c       ppos       : real*8  : position and data (dim (nd+1,*))      
c       xyzf       : real*8  : fit to data (dim nl)
c  output:
c       diff       : real*8  : two element array
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_write_all_fit_data(fname, nd, np,
     >                                    imin, nmin, nmax,
     >                                    ppos, xyzf, diff)
c
        implicit none
c
        include 'mpif.h'
c
        character fname*100
        integer nd, np, imin, nmin, nmax
        real*8 ppos(nd+1, *), xyzf(*), diff(*)        
c
        integer i, j, rank, ierr
        integer fhandle, nwritten, nlocreals
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
     >                       MPI_MODE_CREATE + MPI_MODE_WRONLY,
     >                       MPI_INFO_NULL,
     >                       fhandle, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            file_open = .true.
          else
            mpifunc = 'MPI_File_open'
            exit
          endif
c        
          disp = (imin-1)*FIT_VALUES_PER_PT*SIZE_OF_REAL
c
          call MPI_File_set_view(fhandle, disp,
     >                           MPI_DOUBLE_PRECISION,
     >                           MPI_DOUBLE_PRECISION,
     >                           'native', 
     >                           MPI_INFO_NULL, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_set_view'
            exit 
          endif
c
          nlocreals = (nmax - nmin + 1)*FIT_VALUES_PER_PT
          allocate(buf(nlocreals))
          buf_allocated = .true.
c
          j = 0
          diff(1) = 0.0d0
          diff(2) = 0.0d0
          do i = nmin, nmax
            buf(j+1:j+8) = ppos(1:8,i)
            buf(j+9)     = xyzf(i)
            buf(j+10)    = buf(j+8) - buf(j+9)
            
            diff(1) = diff(1) + buf(j+10)
            diff(2) = diff(2) + (buf(j+10))**2

            j = j+FIT_VALUES_PER_PT
          enddo
c
          call MPI_File_write_all(fhandle, buf, nlocreals,
     >                            MPI_DOUBLE_PRECISION,
     >                            fstat, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_write_all'
            exit
          endif
c
          call MPI_Get_count(fstat, MPI_DOUBLE_PRECISION, 
     >                       nwritten, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            if (nwritten .ne. nlocreals) then
              if (rank .eq. 0) then
                write(*,*) 'Error, unable to write ', nlocreals,
     >                     ' reals to fit data file.'
              endif
            endif
          else
            mpifunc = 'MPI_Get_count'
            exit
          endif
c
          exit
        enddo
c
        if (ierr .ne. MPI_SUCCESS) then
          if (rank .eq. 0) then
            write(*,*) 'Error ', mpifunc, '() returned ', ierr
          endif
        endif
c
        if (buf_allocated) then
          deallocate(buf)
        endif
c
        if (file_open) then
          call MPI_File_close(fhandle, ierr)
        endif
c
        end subroutine mpi_write_all_fit_data
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
c    subroutine sub_sph_wmam_l
c
c       V. Lesur  16/09/2006
c
c   That is for a linearized inversion of the lithosphere field.
c   This subroutine is called from either sub_base_i() or fun_base_f().
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine sub_sph_wmam_l(nub, nb, bc, bp, be, bedotbc, fun)
c
       implicit none
c
       integer nub, nb
       real*8 bc(*), bp(*), be(*)
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
       real*8 bedotbc 
       integer i 
       logical fun 

c 
c    if 'fun' is set to .true. then this routine is being called from fun_base_f() 
c    and returns the dot product of 'be' and 'bc' in 'bedotbc', but 'be' does not
c    contain valid values
c
c    if 'fun' is set to .false. then this routine is being called from sub_base_i()
c    and returns the array 'be', and 'bedotbc' does not contain a valid value.   
c     
c    calculate internal field component
       if (nub .gt. nldp) then
c    this is a sampling point
         call XYZsph_bi0_sample(ilg, rag, bp, d2a, bc,
     >                          dx, dy, dz)

         bex = bp(5)
         bey = bp(6)
         bez = bp(7)
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

       else
c this is a data point 

         bex = bp(5)
         bey = bp(6)
         bez = bp(7)

       endif 

       if (fun) then 
c        computes 'bedotbc' as dot product of 'be' and 'bc' 
         call XYZsph_bi0_fun(ilg, rag, bp, d2a, bc,
     >                       bedotbc, bex, bey, bez)
       else
c        computes 'be' 
         call XYZsph_bi0_sub(ilg, rag, bp, d2a,
     >                       be, bex, bey, bez)
       end if 
c
       return
       end subroutine sub_sph_wmam_l
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sub_base_i
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       subroutine sub_base_i(nub, nb, bc, bp, be)
c
       implicit none
c
       integer nub, nb
       real*8 bc(*), bp(*), be(*)
       real*8 dummybedotbc
c
       call sub_sph_wmam_l(nub, nb, bc, bp, be, dummybedotbc, .false.)
c 
       return
       end subroutine sub_base_i
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    function fun_base_f
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       real*8 function fun_base_f(nub, nb, bc, bp)
c
       implicit none
c
       integer nub, nb
       real*8 bc(*), bp(*), bedotbc
       real*8 dummybe(1) 
c
       call sub_sph_wmam_l(nub, nb, bc, bp, dummybe, bedotbc, .true.)
c
       fun_base_f = bedotbc 
c
       return
       end function fun_base_f
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
