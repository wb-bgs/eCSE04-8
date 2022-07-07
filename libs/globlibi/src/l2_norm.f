ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Function l2_norm
c                  V. Lesur  
c
c       14.07.2011  nub added to parameter list
c
c       Build from the misfit vector MV compute the misfit function
c       for the probability density function associated the L2-norm 
c
c     input:
c       nub     type of data
c	mv	misfit vector
c       i	Element of the misfit vector to use
c    
c     output:  misfit function value
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 function l2_norm(i,nub,mv) 
c
        implicit none
c
        integer i,nub(*)
        real*8 mv(*)
c
        l2_norm=1.d0
c
        return
        end
