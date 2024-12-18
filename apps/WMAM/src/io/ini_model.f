ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_ini_model
c
c  Read in starting model.
c  Rank zero reads in file and then broadcasts data to all other MPI ranks.
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
        integer, parameter :: SIZE_OF_REAL = 8
c
        character fname*100
        integer nparams
        real*8 bc(1:nparams)
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
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_all_ini_model
c
c  Read in starting model collectively.
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
        integer, parameter :: SIZE_OF_REAL = 8
c
        character fname*100
        integer nparams
        real*8 bc(1:nparams)
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
c	
          if (nreals .ne. nparams) then
            if (rank .eq. 0) then
              write(*,*) 'Error, model file should have',
     >                   nparams, ' reals.'
            endif
            nparams = 0
            exit
          endif
c		
          array_of_sizes(1) = nparams
          array_of_subsizes(1) = array_of_sizes(1)
          array_of_starts(1) = 0
c		
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
c		
          call MPI_Type_commit(subarray, ierr)
          if (ierr .eq. MPI_SUCCESS) then
            subarray_type_committed = .true.
          else
            mpifunc = 'MPI_Type_commit'
            exit 
          endif
c			  
          call MPI_File_set_view(fhandle, disp,
     >                           MPI_DOUBLE_PRECISION,
     >                           subarray, 'native', 
     >                           MPI_INFO_NULL, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_set_view'
            exit
          endif
c		
          call MPI_File_read_all(fhandle, bc, 1,
     >                           subarray, fstat, ierr)
          if (ierr .ne. MPI_SUCCESS) then
            mpifunc = 'MPI_File_read_all'
            exit
          endif
c		
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
c
          exit
c
        enddo
c		
        if (ierr .ne. MPI_SUCCESS) then
          if (rank .eq. 0) then
            write(*,*) 'Error ', mpifunc, '() returned ', ierr
          endif
          nparams = 0
        endif
c		
        if (subarray_type_committed) then
          call MPI_Type_free(subarray, ierr)
        endif
c		
        if (file_open) then
          call MPI_File_close(fhandle, ierr)
        endif
c
        end subroutine mpi_read_all_ini_model