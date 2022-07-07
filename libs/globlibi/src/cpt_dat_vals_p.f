cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals_p
c		Vincent Lesur 09/02/2005
c
c       Modified 12.09.2007 V.Lesur 
c       replace common by MPI_Comm_size and MPI_Comm_rank
c       proc_ip(i)=ip changed to proc_ip(i)=ip+1, + adjustment
c
c       Modified to replace the old segmentation scheme by the
c       subroutine thread_segmenter.f (lesur 13.04.2010)
c
c       Parallel interface for cpt_dat_vals.f
c
c       Called: cpt_dat_vals, MPI_ALLGATHERV
c
c       input:
c         nd            Space dimension
c         npm           maximum number of data points
c         np            number of data points
c         ppos(nd+1,*)  point position in ndD
c         nt(*)         basis number for each point
c         nb            Number of base functions
c         BC            base coefficients
c         BS            Base Subroutine to use
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,BC,BS,xyzf)
c
        implicit none
c
        include 'mpif.h'
c
        integer :: nd,npm,npt,nb,nt(*)
        integer :: i,ip,np
        integer, allocatable :: proc_np(:),proc_ip(:)
        real*8  :: ppos(nd+1,*),bc(*),XYZf(*)
c
        external BS
c
        integer ierr,rank,size
        call MPI_Comm_size(MPI_COMM_WORLD,size,ierr)
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
        allocate(proc_np(1:size),proc_ip(1:size))
c
c  All: define what is going to each SP 
        call thread_segmenter(size,npt,proc_np,proc_ip)
c
c  All: find out what are the data to work with
        np=proc_np(rank+1)
        ip=proc_ip(rank+1)
c
c  All: Now does the work
        call cpt_dat_vals(nd,npm,np,nt(ip),ppos(1,ip),nb,BC,BS,xyzf(ip))

c
c  All: loop over each process to broadcast their results
c       MPI_allgatherv would be much better be require too much memory
        do i=1,size
          np=proc_np(i)
          ip=proc_ip(i)
          call MPI_Bcast(xyzf(ip),np,MPI_DOUBLE_PRECISION,i-1
     >                                             ,MPI_COMM_WORLD,ierr)
        enddo
c
        deallocate(proc_np,proc_ip)
        return
        end
