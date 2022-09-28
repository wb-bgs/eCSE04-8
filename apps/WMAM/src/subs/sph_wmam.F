cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     module sph_wmam
c                       V. Lesur 16/09/2006
c                       M. Bareford 18/08/2022
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
c       ilg : integer  : internal SH degree value (see below)
c       ryg : real*8   : reference year for the model
c       rag : real*8   : radius for lithosphere
c
c       dw_fname       : string  : name of burst buffer file holding
c                                  dw values for all data points
c       dw_fhandle     : integer : MPI file handle
c       dw_rank        : integer : MPI rank
c       dw_file_open   : logical : true if MPI dw data file is open
c       dw_offset      : integer : rank offset for dw burst buffer file
c       dw_nlocpts     : integer : number of points assigned to rank
c       dw_nlocdatpts  : integer : number of observed points assigned to rank
c       dw_nlocptsmax  : integer : maximum number of points assigned to rank
c       dw_buf_type_committed : logical : flag indicating if dw_buf_type has been committed
c       dw_buf_type    : integer : MPI dw_buf type
c       dw_buf_size    : integer : number of elements (reals) in dw_buf
c       dw_buf_nbytes  : integer : number of bytes in dw_buf
c       dw_file_nbytes : integer : number of bytes in dw data file
c       dw_buf         : real*8  : allocatable array holding dw values
c                                  for one point
c       dw_offset      : integer : rank offset in units of dw_buf_type
c                                  for dw burst buffer file
c
c    Subroutines:
c       init_sph_wmam()      initialise the global variables
c       precalc_sph_wmam()   precalc dw values and store in burst buffer
c       prepare_dw_read()    open dw data file for reading
c       fini_sph_wmam()      deallocate dw array
c       mpi_read_ref_model()
c       mpi_read_model()
c       mpi_read_data()
c       mpi_read_dw_data()
c       mpi_write_dw_data()
c       sub_sph_wmam()       model subroutine
c       sub_sph_wmam_l()     linearized model subroutine
c       wmam_norm()
c       wmam_var
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        module sph_wmam
c
        implicit none
c
        integer, parameter :: NUM_OF_PTCOMPS = 4
        integer, parameter :: NUM_OF_DWCOMPS = 3
        integer, parameter :: SIZE_OF_REAL = 8
        integer, parameter :: SIZE_OF_PT = NUM_OF_PTCOMPS*SIZE_OF_REAL
c        
        character(len=100) :: dw_fpath
        character(len=100) :: dw_fname
c
        integer ilg
        integer dw_fhandle, dw_rank
        integer dw_nlocdatpts, dw_nlocpts, dw_nlocptsmax
        integer dw_offset
        logical dw_file_open
        logical dw_buf_type_committed
        integer dw_buf_type, dw_buf_size, dw_buf_nbytes
        integer(kind=8) :: dw_file_nbytes
        integer dwx_imin, dwx_imax
        integer dwy_imin, dwy_imax
        integer dwz_imin, dwz_imax
c
        real*8, allocatable :: dw_buf(:) 
c
        real*8 ryg, rag
c
        contains
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine init_sph_wmam
c
c       V. Lesur 16/09/2006
c       M. Bareford 18/08/2022
c
c  Set ilg, ryg (based on 1990) and rag.
c  Set the dw globals.
c
c  input:
c       shdeg       : integer : spherical harmonic degree
c       nparams     : integer : number of parameters
c       nlocdatpts  : integer : number of observed points assigned to rank
c       nlocpts     : integer : number of points assigned to rank
c       nlocptsmax  : integer : maximum number of points assigned to rank
c       npts        : integer : total number of points
c       imin_locpts : integer : first index in to global point array for rank
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine init_sph_wmam(shdeg, nparams, nlocdatpts,
     >                           nlocpts, nlocptsmax,
     >                           npts, imin_locpts)
c
        implicit none
c
        include 'mpif.h'
c
        integer shdeg, nparams
        integer nlocpts, nlocdatpts, nlocptsmax
        integer npts, imin_locpts
c
        character(len=32) :: username
        character(len=5) :: shdeg_str
        integer ierr
c
        call MPI_Comm_rank(MPI_COMM_WORLD, dw_rank, ierr)
c
        call dy2mjd(1990.d0,ryg)
        rag=6371.2d0
        ilg=shdeg
c
        write (shdeg_str, '(I4.4)') shdeg
c
        call GET_ENVIRONMENT_VARIABLE("USER", username)
        dw_fpath='/mnt/lustre/a2fs-nvme/'//TRIM(username)//
     >           '/'//TRIM(shdeg_str)
        dw_fname=TRIM(dw_fpath)//'/dw.dat.bin'
c
        dw_nlocdatpts=nlocdatpts
        dw_nlocpts=nlocpts
        dw_nlocptsmax=nlocptsmax
        dw_offset=imin_locpts-1
        dw_buf_size=NUM_OF_DWCOMPS*nparams
        dw_buf_nbytes=dw_buf_size*SIZE_OF_REAL
        dw_file_nbytes=npts
        dw_file_nbytes=dw_file_nbytes*dw_buf_nbytes
c
        dwx_imin = 1
        dwx_imax = nparams 
        dwy_imin = dwx_imax+1
        dwy_imax = dwx_imax+nparams
        dwz_imin = dwy_imax+1
        dwz_imax = dwy_imax+nparams

        if (dw_rank .eq. 0) then
          write(*,*) ''
          write(*,*) ''
          write(*,*) 'Initialising dw variables'
          write(*,*) 'dw_fname: ', dw_fname
          write(*,*) 'dw_nlocdatpts: ', dw_nlocdatpts
          write(*,*) 'dw_nlocpts: ', dw_nlocpts
          write(*,*) 'dw_nlocptsmax: ', dw_nlocptsmax
          write(*,*) 'dw_offset: ', dw_offset
          write(*,*) 'dw_buf_size: ', dw_buf_size
          write(*,*) 'dw_buf_nbytes: ', dw_buf_nbytes
          write(*,*) 'dw_file_nbytes: ', dw_file_nbytes
          write(*,*) ''
          write(*,*) 'dwx_imin: ', dwx_imin
          write(*,*) 'dwx_imax: ', dwx_imax
          write(*,*) 'dwy_imin: ', dwy_imin
          write(*,*) 'dwy_imax: ', dwy_imax
          write(*,*) 'dwz_imin: ', dwz_imin
          write(*,*) 'dwz_imax: ', dwz_imax
          write(*,*) ''
          write(*,*) ''
        endif
c
        dw_buf_type_committed = .false.
        dw_file_open = .false.
c
        call MPI_Type_contiguous(dw_buf_size,
     >                           MPI_DOUBLE_PRECISION,
     >                           dw_buf_type, ierr)
        if (ierr .ne. MPI_SUCCESS) then
          write(*,*) dw_rank,
     >               ': Error MPI_Type_contiguous() returned ',
     >               ierr
          stop
        endif
c
        call MPI_Type_commit(dw_buf_type, ierr)
        if (ierr .eq. MPI_SUCCESS) then
          dw_buf_type_committed = .true.
        else
          write(*,*) dw_rank,
     >               ': Error MPI_Type_commit() returned ',
     >               ierr
          stop
        endif
c        
        if (dw_buf_type_committed) then
          allocate(dw_buf(dw_buf_size))
        endif
c
        end subroutine init_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine precalc_sph_wmam
c
c  Calculate and write dw values to burst buffer.
c
c  input:
c       nd        : integer : nd+1 is the lead dim of ppos
c       nparams   : integer : number of parameters
c       ppos(*,*) : real*8  : array of data values
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine precalc_sph_wmam(nd, nparams,  ppos)
c
        implicit none
c
#ifdef CRAYPAT
        include 'pat_apif.h'
#endif
c
        include 'mpif.h'
c
        integer nd, nparams
        real*8 ppos(nd+1,*)
c
        logical fexists
        integer ierr, i
        integer(kind=8) :: fsize = 0
        character(len=32) :: fsize_str
        character(len=200) :: cmd  
c
#ifdef CRAYPAT
        integer cpstat
#endif
c
        dw_file_open = .false.
c
        inquire(FILE=dw_fname, EXIST=fexists)
c
        if (fexists) then
          inquire(FILE=dw_fname, SIZE=fsize)
          if (fsize .eq. dw_file_nbytes) then
c  correctly-sized dw data file already exists
c  assume pre-calculations have already been written
            if (dw_rank .eq. 0) then
              write(*,*) 'Found existing dw data file at ',
     >                   dw_fname
              write(*,*) ''
            endif
            return
          else
c  delete incorrectly sized dw data file
            cmd='rm -rf '//dw_fpath
            call EXECUTE_COMMAND_LINE(cmd)
          endif
        endif
c
c  create file path to dw data file
        cmd='mkdir -p '//TRIM(dw_fpath)
        call EXECUTE_COMMAND_LINE(cmd)
c
c  set the strip count to -1 (i.e., use all OSTs)
        cmd='lfs setstripe -c -1 '//TRIM(dw_fpath)
        call EXECUTE_COMMAND_LINE(cmd)
c
c  set the size of the dw data file
        fsize = dw_file_nbytes
        write(fsize_str,*) fsize
        cmd='truncate -s '//TRIM(fsize_str)//' '//TRIM(dw_fname)
        call EXECUTE_COMMAND_LINE(cmd)
c
c  open dw data file for writing
        call MPI_File_open(MPI_COMM_WORLD, dw_fname,
     >                     MPI_MODE_WRONLY, MPI_INFO_NULL,
     >                     dw_fhandle, ierr)
        if (ierr .eq. MPI_SUCCESS) then
          dw_file_open = .true.
        else
          write(*,*) dw_rank,
     >               ': Error MPI_File_open(WRONLY) returned ',
     >               ierr
          stop
        endif
c 
c  data points
        do i=1,dw_nlocpts
          dw_buf(dwx_imin:dwz_imax)=0.0d0
c
#ifdef CRAYPAT
          call PAT_record(PAT_STATE_ON,cpstat)
          call PAT_region_begin(1,"XYZsph_bi0",cpstat)
#endif
c
c  calculate internal field component
          call XYZsph_bi0(ilg,rag,ppos(1,i),
     >                    dw_buf(dwx_imin:dwx_imax),
     >                    dw_buf(dwy_imin:dwy_imax),
     >                    dw_buf(dwz_imin:dwz_imax))
c
#ifdef CRAYPAT
          call PAT_region_end(1,cpstat)
          call PAT_record(PAT_STATE_OFF,cpstat)
#endif
c
c  write components to burst buffer 
          call mpi_write_dw_data(i, .false.)
        enddo
c
c  do dummy writes so that mpi collective write completes
        do i=dw_nlocpts+1,dw_nlocptsmax
          call mpi_write_dw_data(dw_nlocpts, .false.)
        enddo
c
c  close dw data file after writing
        if (dw_file_open) then
          call MPI_File_close(dw_fhandle, ierr)
          dw_file_open = .false.
        endif
c
        end subroutine precalc_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine prepare_dw_read
c
c  Open dw data file for reading.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine prepare_dw_read()
c
        implicit none
c
        include 'mpif.h'
c
        integer ierr
c
        dw_file_open = .false.
c
        call MPI_File_open(MPI_COMM_WORLD, dw_fname,
     >                     MPI_MODE_RDONLY, MPI_INFO_NULL,
     >                     dw_fhandle, ierr)
        if (ierr .eq. MPI_SUCCESS) then
          dw_file_open = .true.
        else
          write(*,*) dw_rank,
     >               ': Error MPI_File_open(RDONLY) returned ',
     >               ierr
          stop
        endif
c
        end subroutine prepare_dw_read
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine fini_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine fini_sph_wmam()
c
        implicit none
c
        include 'mpif.h'
c
        integer ierr
c
c  close dw data file after reading
        if (dw_file_open) then
          call MPI_File_close(dw_fhandle, ierr)
          dw_file_open = .false.
        endif
c
        deallocate(dw_buf)
c
        if (dw_buf_type_committed) then
          call MPI_Type_free(dw_buf_type, ierr)
          dw_buf_type_committed = .false.
        endif
c
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

        end
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

        end
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
c
          nlocreals = nlocdatpts*NUM_OF_PTCOMPS
          allocate(buf(nlocreals))
          buf_allocated = .true.

          disp = (imin_locdatpts-1)*SIZE_OF_PT
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
          call MPI_File_read_all(fhandle, buf, nlocreals,
     >                           MPI_DOUBLE_PRECISION,
     >                           fstat, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_read_all'
            exit
          endif
c
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

        end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_dw_data
c
c  Read in dw values associated with a specific point.
c
c  input:
c       ipt : integer : one-based index identifying a particular
c                       point assigned to this rank
c       debug : logical : if true check read count and if any point
c                         read is nan
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_dw_data(ipt, debug)
c
        implicit none
c
        include 'mpif.h'
c
	integer ipt
        logical debug
c
        integer ierr, i
        integer :: nread = 0
        integer (kind=MPI_OFFSET_KIND) :: disp = 0
        integer fstat(MPI_STATUS_SIZE)
c
        disp = dw_offset + (ipt-1)
        disp = disp*dw_buf_nbytes
c
        call MPI_File_set_view(dw_fhandle, disp,
     >                         MPI_DOUBLE_PRECISION,
     >                         dw_buf_type, 'native', 
     >                         MPI_INFO_NULL, ierr)
        if (ierr .ne. MPI_SUCCESS) then
          write(*,*) dw_rank,
     >               ': Error, MPI_File_set_view() returned ',
     >               ierr
          stop 
        endif
c
        call MPI_File_read_all(dw_fhandle, dw_buf, 1,
     >                         dw_buf_type, fstat, ierr)
        if (ierr .ne. MPI_SUCCESS) then
          write(*,*) dw_rank,
     >               ': Error, MPI_File_read_all() returned ',
     >               ierr
          stop
        endif
c
        if (debug) then
          call MPI_Get_count(fstat, dw_buf_type, nread, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            if (nread .ne. 1) then
              write(*,*) 'Error, unable to read ',
     >                    dw_buf_size, ' reals from dw file'
              stop
            endif
          else
            write(*,*) dw_rank,
     >                 ': Error, MPI_Get_count() returned ',
     >                 ierr
            stop
          endif
c
          do i=1,dw_buf_size
            if (ISNAN(dw_buf(i))) then
              write(*,*) dw_rank,
     >                   ': Read error, dw value ', i,
     >                   ' for point ', ipt, ' is nan'
              stop 
            endif
          enddo
        endif
c
        end subroutine mpi_read_dw_data
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_write_dw_data
c
c  Write the dw values associated with a specific point.
c
c  input:
c       ipt : integer : one-based index identifying a particular
c                       point assigned to this rank
c       debug : logical : if true check write count is correct
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_write_dw_data(ipt, debug)
c
        implicit none
c
        include 'mpif.h'
c
        integer ipt
        logical debug
c
        integer ierr
        integer :: nwritten = 0
        integer (kind=MPI_OFFSET_KIND) :: disp = 0
        integer fstat(MPI_STATUS_SIZE)
c        
        disp = dw_offset + (ipt-1)
        disp = disp*dw_buf_nbytes
c
        call MPI_File_set_view(dw_fhandle, disp,
     >                         MPI_DOUBLE_PRECISION,
     >                         dw_buf_type, 'native', 
     >                         MPI_INFO_NULL, ierr)
        if (ierr .ne. MPI_SUCCESS) then
          write(*,*) dw_rank,
     >               ': Error, MPI_File_set_view() returned ',
     >               ierr
          stop 
        endif
c
        call MPI_File_write_all(dw_fhandle, dw_buf, 1,
     >                          dw_buf_type, fstat, ierr)
        if (ierr .ne. MPI_SUCCESS) then
          write(*,*) dw_rank,
     >               ': Error, MPI_File_write_all() returned ',
     >               ierr
          stop
        endif
c
        if (debug) then
          call MPI_Get_count(fstat, dw_buf_type, nwritten, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            if (nwritten .ne. 1) then
              write(*,*) dw_rank, ': Error, unable to write ',
     >                   dw_buf_size, ' reals to dw file'
              stop
            endif
          else
            write(*,*) dw_rank,
     >                 ': Error MPI_Get_count() returned ',
     >                 ierr
            stop
          endif
        endif
c
        end subroutine mpi_write_dw_data
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sub_sph_wmam
c
c       V. Lesur  16/09/2006
c
c   That is for a full non-linear inversion of the lithosphere field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sub_sph_wmam(iof, nub, nb, bc, bp, be)
c
        implicit none
c
        character iof
        integer nub, nb
        real*8 bc(*), bp(*), be(*)
c
        integer i
        real*8 dx,dy,dz,dd
        real*8 xy_c, xz_c, yx_c, yz_c, zx_c, zy_c
        real*8 bex, bey, bez
        real*8 bex2, bey2, bez2
c
c
        be(1:nb) = 0.0d0

c  calculate internal field component
c  read in dw for point nub from burst buffer
        call mpi_read_dw_data(nub, .false.)
c
        if (nub .gt. dw_nlocdatpts) then
c
c  sampling point
          dx = dot_product(dw_buf(dwx_imin:dwx_imax), bc(1:nb))
          dy = dot_product(dw_buf(dwy_imin:dwy_imax), bc(1:nb))
          dz = dot_product(dw_buf(dwz_imin:dwz_imax), bc(1:nb))
c
c  calculate field unit vector (lithos)
          dd = dsqrt(bp(5)**2 + bp(6)**2 + bp(7)**2)
          bex = bp(5) / dd
          bey = bp(6) / dd
          bez = bp(7) / dd
c
          xy_c = dx*bey - bex*dy
          xz_c = dx*bez - bex*dz
          yx_c = dy*bex - bey*dx
          yz_c = dy*bez - bey*dz
          zx_c = dz*bex - bez*dx
          zy_c = dz*bey - bez*dy
c
          dd = dsqrt(yz_c**2 + xz_c**2 +xy_c**2)
c
          if (iof.ne.'f') then
            bex2 = (xz_c*bez + xy_c*bey) / dd
            bey2 = (yz_c*bez + yx_c*bex) / dd
            bez2 = (zy_c*bey + zx_c*bex) / dd
c
            be(1:nb) = dw_buf(dwx_imin:dwx_imax)*bex2
     >               + dw_buf(dwy_imin:dwy_imax)*bey2
     >               + dw_buf(dwz_imin:dwz_imax)*bez2
          else
c  full non-linear forward
            be(1:1) = dd
          endif
        
        else
c  regular data
c
c  calculate full field (lithos + core)
          dx = bp(5) + dot_product(dw_buf(dwx_imin:dwx_imax), bc(1:nb))
          dy = bp(6) + dot_product(dw_buf(dwy_imin:dwy_imax), bc(1:nb))
          dz = bp(7) + dot_product(dw_buf(dwz_imin:dwz_imax), bc(1:nb))

          dd = dsqrt(dx**2 + dy**2 + dz**2)
c
          if (iof.ne.'f') then
            bex = dx / dd
            bey = dy / dd
            bez = dz / dd
c
            be(1:nb) = dw_buf(dwx_imin:dwx_imax)*bex
     >               + dw_buf(dwy_imin:dwy_imax)*bey
     >               + dw_buf(dwz_imin:dwz_imax)*bez
          else
c  full non-linear forward
            dx = dsqrt(bp(5)**2 + bp(6)**2 + bp(7)**2)
            be(1:1) = dd - dx
          endif
        
        endif
c
c  do dummy reads so that collective mpi reads complete
        if (nub .eq. dw_nlocpts) then
          do i=dw_nlocpts+1,dw_nlocptsmax
            call mpi_read_dw_data(dw_nlocpts, .false.)
          enddo
        endif
c
        end subroutine  sub_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sub_sph_wmam_l
c
c       V. Lesur  16/09/2006
c
c   That is for a linearized inversion of the lithosphere field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sub_sph_wmam_l(iof, nub, nb, bc, bp, be)
c
        implicit none
c
        character iof
        integer nub, nb
        real*8 bp(*), be(*), bc(*)
c
        integer i
        real*8 dx, dy, dz, dd
        real*8 xy_c, xz_c, yx_c, yz_c, zx_c, zy_c
        real*8 bex, bey, bez
        real*8 bex2, bey2, bez2
c
        be(1:nb) = 0.0d0

c    calculate full field (core)
        dd = dsqrt(bp(5)**2 + bp(6)**2 + bp(7)**2)
        bex = bp(5) / dd
        bey = bp(6) / dd
        bez = bp(7) / dd

c    calculate internal field component
c    read in dw for point nub from burst buffer
        call mpi_read_dw_data(nub, .false.)
c
        if (nub .gt. dw_nlocdatpts) then
c    sampling point
          dx = dot_product(dw_buf(dwx_imin:dwx_imax), bc(1:nb))
          dy = dot_product(dw_buf(dwy_imin:dwy_imax), bc(1:nb))
          dz = dot_product(dw_buf(dwz_imin:dwz_imax), bc(1:nb))
c
          xy_c = dx*bey - bex*dy
          xz_c = dx*bez - bex*dz
          yx_c = dy*bex - bey*dx
          yz_c = dy*bez - bey*dz
          zx_c = dz*bex - bez*dx
          zy_c = dz*bey - bez*dy
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

        dw_buf(dwx_imin:dwx_imax) = dw_buf(dwx_imin:dwx_imax) * bex
        dw_buf(dwy_imin:dwy_imax) = dw_buf(dwy_imin:dwy_imax) * bey
        dw_buf(dwz_imin:dwz_imax) = dw_buf(dwz_imin:dwz_imax) * bez

        be(1:nb) = dw_buf(dwx_imin:dwx_imax)
     >           + dw_buf(dwy_imin:dwy_imax)
     >           + dw_buf(dwz_imin:dwz_imax)

        if (iof .eq. 'f') then
          be(1:nb) = be(1:nb) * bc(1:nb)
        endif
c
c  do dummy reads so that collective mpi reads complete
        if (nub .eq. dw_nlocpts) then
          do i=dw_nlocpts+1,dw_nlocptsmax
            call mpi_read_dw_data(dw_nlocpts, .false.)
          enddo
        endif
c
        end subroutine  sub_sph_wmam_l
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
c          wmam_norm=lesur_norm(dc,ae,i,mv)
         case default
           wmam_norm=l2_norm(i,nub,mv)
c          wmam_norm=l1_norm(i,nub,mv)
         end select
c
         return
         end function wmam_norm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c              vv=lesur_var2(dc,ae,j,k,dv,mv(i),wgh(i))
             case default
               vv=l2_var(nub(i),j,k,dv,mv(i),wgh(i))
c              vv=l1_var(nub(i),j,k,dv,mv(i),wgh(i))
             end select
c
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
