cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals_p2
c		Vincent Lesur 09/02/2005
c
c       Parallel interface for cpt_dat_vals.f
c
c       Called: cpt_dat_vals_p2
c
c       input:
c         nd            Space dimension
c         nlocpts       number of data+sampling points local to rank
c         ppos(nd+1,*)  point position in ndD
c         nb            Number of base functions
c         bc            base coefficients
c         fun_base      Base function to use
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p2(nd, nlocpts, ppos,
     >                             nb, bc, fun_base, xyzf)
c
        implicit none
c
        integer :: nd, nb, nlocpts
        real*8  :: ppos(nd+1,*), bc(*)
        real*8  :: fun_base, xyzf(*)
c
        integer :: i
c
        external fun_base

        do i=1,nlocpts
          xyzf(i) = fun_base(i,nb,bc,ppos(1,i))
        enddo

        return
        end