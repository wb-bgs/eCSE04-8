ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_ref_model
c
c  Read in reference model.
c  Rank 0 does the read and then broadcasts data to all other ranks.
c
c  input:
c       fname   : character : name of reference model file
c       ncoeffs : integer : number of coefficients/parameters
c  output:
c       ryg     : real*8  : reference year for the model
c       bc(*)   : real*8  : coefficients/parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_ref_model(fname, ncoeffs, ryg, bc)
c
        implicit none
c
        include 'mpif.h'
c
        integer, parameter :: SIZE_OF_REAL = 8
c
        character fname*100
        integer ncoeffs
        real*8 ryg
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
c  globlibi subroutine
        external dy2mjd
c
c
        nreals = 1 + ncoeffs
        allocate(buf(1:nreals))
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_all_ref_model
c
c  Read in reference model collectively.
c  All ranks read entire file via MPI_File_read_all().
c
c  input:
c       fname   : character : name of reference model file
c  output:
c       ncoeffs : integer : number of coefficients/parameters
c       ryg     : real*8  : reference year for the model
c       bc(*)   : real*8  : coefficients/parameters
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_all_ref_model(fname, ncoeffs, ryg, bc)
c
        implicit none
c
        include 'mpif.h'
c
        integer, parameter :: SIZE_OF_REAL = 8
c
        character fname*100
        integer ncoeffs
        real*8 ryg
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
c  globlibi subroutine
        external dy2mjd
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

          allocate(buf(1:nreals))
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
