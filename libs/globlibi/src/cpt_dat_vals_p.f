cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals_p
c		Vincent Lesur 09/02/2005
c
c       Parallel interface for cpt_dat_vals.f
c
c       Called: cpt_dat_vals_p
c
c       input:
c         nd            Space dimension
c         np            number of data points
c         nt            basis number for each point
c         ppos(nd+1,*)  point position in ndD
c         nb            Number of base functions
c         bc            base coefficients
c         sub_base      Base Subroutine to use
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p(nd, np, nt, ppos, nb, bc,
     >                            sub_base, xyzf)
c
        implicit none
c
        integer :: nd,np,nt,nb
        real*8  :: ppos(nd+1,*),bc(*),xyzf(*)
c
        integer :: i
        real*8, allocatable :: dw(:)
c
        external sub_base
c
        allocate(dw(1:nb))
c
        do i=1,np
          call sub_base('f',nt,nb,bc,ppos(1,i),dw)
          xyzf(i) = SUM(dw)
        enddo
c
        deallocate(dw)
c
        return
        end