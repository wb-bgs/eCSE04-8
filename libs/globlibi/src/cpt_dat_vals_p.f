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
c         proc_np(*)    block lengths
c         proc_ip(*)    data pointers
c         nt(*)         basis number for each point
c         ppos(nd+1,*)  point position in ndD
c         nb            Number of base functions
c         bc            base coefficients
c         sub_base      Base Subroutine to use
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p(nd, proc_np, proc_ip, nt,
     >                            ppos, nb, bc, sub_base, xyzf)
c
        implicit none
c
        include 'mpif.h'
c
        integer :: nd,nb,nt(*)
        integer :: ip,np
        integer :: proc_np(*),proc_ip(*)
        real*8  :: ppos(nd+1,*),bc(*),XYZf(*)
c
        external sub_base
c
        integer ierr, rank
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
c  All: find out what are the data to work with
        np = proc_np(rank+1)
        ip = proc_ip(rank+1)+1
c
c  All: Now does the work
        call cpt_dat_vals(nd, np, nt(ip), ppos(1,ip), nb,
     >                    bc, sub_base, xyzf(ip))
 
        call MPI_ALLGATHERV(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
     >                      xyzf, proc_np, proc_ip,
     >                      MPI_DOUBLE_PRECISION,
     >                      MPI_COMM_WORLD, ierr)

        return
        end