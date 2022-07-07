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
c         npm           maximum number of data points
c         np            number of data points
c	  ppos(nd+1,*)	point position in ndD
c         nt(*)         basis number for each point
c         nb		Number of base functions
c         BC            base coefficients
c         BS            Base Subroutine to use
c
c	output:
c         XYZF(*)	X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine cpt_dat_vals(nd,npm,np,nt,ppos,nb,BC,BS,xyzf)
c
        implicit none
c
        integer :: i,nd,npm,np,nb,nt(*)
        real*8 :: ppos(nd+1,*),bc(*),XYZf(*)
        real*8, allocatable :: dw(:)
c
        external BS
c
        allocate(dw(1:nb))
        do i=1,np
          call cpt_dat(nt(i),nb,BS,BC,ppos(1,i),xyzf(i),dw)
        enddo
        deallocate(dw)
c
        return
        end
