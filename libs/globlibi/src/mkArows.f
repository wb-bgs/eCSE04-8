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
c         issampt       true of first point is sampling point 
c         shdeg         max SH degree value
c         nb            number of base functions
c         nd		space dimension
c         np		number of data points
c         npmax         number max of data points handled together
c         d2a           pre-computed array for mklf_F2()
c         bc            estimation of base function coefficients
c         ppos		point position in ndD
c		
c       output:
c	  aa(NP,NB)	matrix of conditions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine mkArows(issampt, shdeg, nb, nd, np, npmax,
     >                     d2a, bc, ppos, aa)
c
        implicit none
c
        logical issampt
	integer shdeg, nb, nd, np, npmax
        real*8 d2a(0:shdeg), bc(nb), ppos(nd+1,np)
        real*8 aa(nb,npmax)
c
        integer i
c
        external wmam_sub
!$omp declare target
c
c
        do i=1,np
          call wmam_sub(issampt,
     >                  shdeg, nb, nd,
     >                  d2a, bc, ppos(1,i),
     >                  aa(1,i))
        enddo
c
c
        return
	end