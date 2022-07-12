cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals
c		Vincent Lesur 09/02/2005
c
c       calling recusive subroutine, should be compiled using F90
c       previously called mkxyzf.
c
c       for a given set of basis coefficients BC and basis functions
c       calculate the mainfield X,Y,Z or F components (depending on nt value)
c       for a set of data point in a nd dimmensional space
c
c       Called: cpt_dat
c
c       input:
c         nd	        Space dimension
c         np            number of data points
c	  ppos(nd+1,*)	point position in ndD
c         nt(*)         basis number for each point
c         nb		Number of base functions
c         bc            base coefficients
c         sub_base      Base Subroutine to use
c
c	output:
c         xyzf(*)	X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine cpt_dat_vals(nd,np,nt,ppos,nb,bc,sub_base,xyzf)
c
        implicit none
c
        integer :: i,nd,np,nb,nt(*)
        real*8 :: ppos(nd+1,*),bc(*),xyzf(*)
        real*8, allocatable :: dw(:)
c
        external sub_base
c
        allocate(dw(1:nb))
        do i=1,np
          call cpt_dat(nt(i),nb,sub_base,bc,ppos(1,i),xyzf(i),dw)
        enddo
        deallocate(dw)
c
        return
        end