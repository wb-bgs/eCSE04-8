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
        integer, parameter :: SIZE_OF_REAL = 8
        integer, parameter :: FIT_VALUES_PER_PT = 10
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
          allocate(buf(1:nlocreals))
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