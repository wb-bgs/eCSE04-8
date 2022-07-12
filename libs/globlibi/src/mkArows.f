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
c         nt(np)	Type of data (1->X,2->Y,3->Z,4->F ...) dim min (NP) 
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
	subroutine mkArows(np,nt,nd,nb,ppos,sub_base,bc,aa)
c
        implicit none
c
	integer np,nt(*),nd,nb
        real*8 ppos(nd+1,*),aa(np,*),bc(*)
        real*8, allocatable :: row(:)
c
        integer ip,ib
c
        external sub_base
        allocate(row(1:nb))
c
        do ip=1,np
          call sub_base('i',nt(ip),nb,bc,ppos(1,ip),row)
          do ib=1,nb
            aa(ip,ib)=row(ib)
          enddo
        enddo
        deallocate(row)
c
        return
	end