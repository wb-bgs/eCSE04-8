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
c         nlocdatpts    number of data points local to rank
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
        subroutine cpt_dat_vals_p2(nd, nlocdatpts, nlocpts, ppos,
     >                             nb, bc, fun_base, xyzf)
c
        implicit none
c
        integer :: nd, nb
        integer :: nlocdatpts, nlocpts
        real*8  :: ppos(nd+1,*), bc(*)
        real*8  :: fun_base, xyzf(*)
c
        integer :: i
        real*8, allocatable :: dw1(:)
c
        external fun_base
c
        allocate(dw1(nb))

c  data points
        do i=1,nlocdatpts
          xyzf(i) = fun_base(1,nb,bc,ppos(1,i),dw1)
        enddo

c  sampling points
        do i=nlocdatpts+1,nlocpts
          xyzf(i) = fun_base(100,nb,bc,ppos(1,i),dw1)
        enddo
c
        deallocate(dw1)
c
        return
        end
