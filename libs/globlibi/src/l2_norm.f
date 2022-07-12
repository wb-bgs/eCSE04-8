ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Function l2_norm
c                  V. Lesur  
c
c       14.07.2011  nub added to parameter list
c
c       Build from the misfit vector MV compute the misfit function
c       for the probability density function associated the L2-norm 
c
c       Actually, this is a stub function that simply returns 1.d0.
c
c     input:
c    
c     output:  misfit function value
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 function l2_norm() 
c
        implicit none
c
        l2_norm=1.d0
c
        return
        end