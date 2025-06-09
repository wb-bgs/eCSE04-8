	subroutine init_mpi_status_types(nb, bc,
     >                                   inv_stat, src_stat,
     >                                   MPI_INVERSION_STATUS,
     >                                   MPI_SEARCH_STATUS)
c
        use mpi
c
        implicit none
c
c       include 'mpif.h'
        include 'mpi_status_types.h'
c
	integer :: nb
	real*8 :: bc(1:nb)
        type(inversion_status) :: inv_stat
        type(search_status) :: src_stat
        integer :: MPI_INVERSION_STATUS
        integer :: MPI_SEARCH_STATUS 

        integer, parameter :: STAT_NBLOCKS=2
        integer(KIND=MPI_ADDRESS_KIND) :: base

        integer, dimension(STAT_NBLOCKS) :: inv_stat_lens
        integer, parameter, dimension(STAT_NBLOCKS)
     >    :: inv_stat_types = (/ MPI_CHARACTER, 
     >                           MPI_DOUBLE_PRECISION /)
        integer(KIND=MPI_ADDRESS_KIND)
     >    :: inv_stat_disps(STAT_NBLOCKS)

        integer, parameter, dimension(STAT_NBLOCKS)
     >      :: src_stat_lens = (/ 1, 1 /)
        integer, parameter, dimension(STAT_NBLOCKS)
     >      :: src_stat_types = (/ MPI_CHARACTER,
     >                             MPI_DOUBLE_PRECISION /)
        integer(KIND=MPI_ADDRESS_KIND)
     >      :: src_stat_disps(STAT_NBLOCKS)

        integer :: ierr
	
c  Commit MPI_INVERSION_STATUS type
        inv_stat_lens = (/ 5, nb /)

        allocate(inv_stat%bc(nb))
        inv_stat%bc(1:nb)=bc(1:nb)
c       deallocate(bc)

        call MPI_GET_ADDRESS(inv_stat%yon, inv_stat_disps(1), ierr) 
	call MPI_GET_ADDRESS(inv_stat%bc, inv_stat_disps(2), ierr)
        
        base=inv_stat_disps(1)
        inv_stat_disps(1)=inv_stat_disps(1)-base
        inv_stat_disps(2)=inv_stat_disps(2)-base

	call MPI_TYPE_CREATE_STRUCT(stat_nblocks, inv_stat_lens,
     >                              inv_stat_disps, inv_stat_types,
     >                              MPI_INVERSION_STATUS, ierr) 
	call MPI_TYPE_COMMIT(MPI_INVERSION_STATUS, ierr)
	    
c  Commit MPI_SEARCH_STATUS type
	call MPI_GET_ADDRESS(src_stat%yon_ct, src_stat_disps(1), ierr) 
        call MPI_GET_ADDRESS(src_stat%stp, src_stat_disps(2), ierr)

        base=src_stat_disps(1)
        src_stat_disps(1)=src_stat_disps(1)-base
        src_stat_disps(2)=src_stat_disps(2)-base

	call MPI_TYPE_CREATE_STRUCT(stat_nblocks, src_stat_lens,
     >                              src_stat_disps, src_stat_types,
     >                              MPI_SEARCH_STATUS, ierr) 
        call MPI_TYPE_COMMIT(MPI_SEARCH_STATUS, ierr)
	
        end subroutine


        subroutine free_mpi_status_types(nb,bc,inv_stat,
     >                                   MPI_INVERSION_STATUS,
     >                                   MPI_SEARCH_STATUS)
c
        use mpi
c
        implicit none
c
c       include 'mpif.h'
        include 'mpi_status_types.h'
c
        integer :: nb
        real*8 :: bc(*)
        type(inversion_status) :: inv_stat
        integer :: MPI_INVERSION_STATUS
        integer :: MPI_SEARCH_STATUS 
c
        integer :: ierr
	     
        call MPI_TYPE_FREE(MPI_INVERSION_STATUS, ierr)
	call MPI_TYPE_FREE(MPI_SEARCH_STATUS, ierr)

c       allocate(bc(1:nb))
        bc(1:nb)=inv_stat%bc(1:nb)
        deallocate(inv_stat%bc)
	     
        end subroutine
