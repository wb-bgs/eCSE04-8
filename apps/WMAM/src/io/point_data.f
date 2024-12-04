ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine mpi_read_all_data
c
c  Read in point data collectively.
c  Read just those points that have been assigned to the MPI rank.
c
c  input:
c       fname      : character : name of data file
c       nd         : integer : nd+1 is the lead dim of ppos
c       nlocpts    : integer : total number of points assigned to rank
c       ndatpts    : integer : number of data points in file
c       nlocdatpts : integer : number of data points assigned to rank
c       ryg        : real*8  : reference year for the model
c       imin_locdatpts : integer : one-based rank index for global
c                                  array of data points
c  output:
c       ppos(*,*) : real*8 : array of data points assigned to rank
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mpi_read_all_data(fname, nd, nlocpts, ndatpts,
     >                               nlocdatpts, imin_locdatpts,
     >                               ryg, ppos)
c
        implicit none
c
        include 'mpif.h'
c
        integer, parameter :: NUM_OF_PTCOMPS = 4
        integer, parameter :: SIZE_OF_REAL = 8
        integer, parameter :: SIZE_OF_PT = NUM_OF_PTCOMPS*SIZE_OF_REAL
c
        character fname*100
        integer nd, nlocpts, ndatpts, nlocdatpts, imin_locdatpts
        real*8 ryg, ppos(1:nd+1,1:nlocpts)
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
          allocate(buf(1:nlocreals))
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