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
c	  NP		Number of data points
c         npm           not used anymore
c         NT(NP)	Type of data (1->X,2->Y,3->Z,4->F ...) dim min (NP) 
c         ND		space dimension
c         nb            Number of base functions
c         ppos		point position in ndD
c         BS            the "Base functions" subroutine to use
c         XYZF          X,Y,Z component and total field dim min (NP,4)
c		
c       output:
c	  AA(NP,NB)	matrix of conditions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mkArows(np,npm,nt,nd,nb,ppos,bs,bc,AA)
c
        implicit none
c
	integer np,npm,nt(*),nd,nb
        real*8 ppos(nd+1,*),AA(np,*),bc(*)
        real*8, allocatable :: row(:)
c
        integer ip,ib,j
c
        external bs
        allocate(row(1:nb))
c
        do ip=1,np
          call BS('i',nt(ip),nb,bc,ppos(1,ip),row)
          do ib=1,nb
            AA(ip,ib)=row(ib)
          enddo
        enddo
        deallocate(row)
c
        return
	end
