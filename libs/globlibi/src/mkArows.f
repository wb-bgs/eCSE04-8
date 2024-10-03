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
c         shdeg         max SH degree value
c         nb            number of base functions
c         nd		space dimension
c         np		number of data points
c         ip            point index
c         nlocdatpts    number of data points assigned to rank
c         d2a           pre-computed array for mklf_F2()
c         bc            estimation of base function coefficients
c         ppos		point position in ndD
c		
c       output:
c	  aa(NP,NB)	matrix of conditions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
        subroutine mkArows(shdeg, nb, nd, np,
     >                     ip, nlocdatpts,
     >                     d2a, bc, ppos,
     >                     aa)
c
        implicit none
c
        integer shdeg, nb, nd, np
        integer ip, nlocdatpts
        real*8 d2a(0:shdeg), bc(nb), ppos(nd+1,np)
        real*8 aa(nb,np)
c
        integer i
c
        external wmam_sub
c
c
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
c
        do i=1,np
          call wmam_sub((ip>nlocdatpts),
     >                  shdeg, nb, nd,
     >                  d2a, bc, ppos(1,i),
     >                  aa(1,i))
          ip = ip + 1
        enddo
c
        return
	end
