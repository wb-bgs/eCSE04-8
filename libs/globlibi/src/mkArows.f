cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	subroutine mkArows
c		V. Lesur 09/02/2005
c
c       previouslu mkcondmat, should be compile using F90
c
c       calculate the condition matrix for NP data points.
c       call recurcive subroutines
c
c	input:
c	  np		Number of data points
c         nt  	        Type of data (1->X,2->Y,3->Z,4->F ...) 
c         ip            starting number for basis number for each point
c                       (1->X,2->Y,3->Z,4->F ...)
c         nd		space dimension
c         nb            Number of base functions
c         ppos		point position in ndD
c         sub_base      the "Base functions" subroutine to use
c         xyzf          X,Y,Z component and total field dim min (NP,4)
c		
c       output:
c	  aa(NP,NB)	matrix of conditions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mkArows(np,ip,nd,nb,ppos,sub_base,bc,aa)
c
        implicit none
c
	integer np,ip,nd,nb
        real*8 ppos(nd+1,*),bc(*)
        real*8 aa(nb,*)
c
        integer i
c
        external sub_base
c
        do i=1,np
          call sub_base(ip+i-1,nb,bc,
     >                  ppos(1,i),
     >                  aa(1,i))
        enddo
c
        return
	end