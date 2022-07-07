ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Function l2_var
c                  V. Lesur	25.11.2011  
c
c     modified on 21.08.2012 such that var & npt are not modified 
c     during a call (V. Lesur)
c
c      depending on input do one of the two:
c      1-  np.ge.0
c            l2_var =  l2_var 
c               + sum_{i=1}^{np} wght(i)*mv(i)**2
c      3-  np.lt.0
c            l2_var = l2_var of the combined (-np) l2_var
c
c       l2_var is computed with an L2 norm
c
c     input np>0:
c       nub     type of data
c       npt     np total number of point used to compute var
c       var     values of var computed with npt points
c       np      number of point added in this itteration
c       mv      misfit vector dimmin np
c       wgh     vector of 1/variances dimmin np
c     input np<0:
c       np      (-1*) number of var to be added
c       mv      var values  dimmin -np
c       wgh     number of point used to compute vars dimmin -np
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 function l2_var(nub,npt,np,var,mv,wgh) 
c
        implicit none
c
        integer npt,np,i,nub(*),nn
        real*8 mv(*),var,wgh(*),dd
c
        if(np.ge.0)then
          dd=var
          do i=1,np
            dd=dd+wgh(i)*mv(i)**2
          enddo
          l2_var=dd
        else
          dd=0.0d0
          nn=0
          do i=1,-np
            dd=dd+mv(i)
            nn=nn+idnint(wgh(i))
          enddo
          l2_var=dd
        endif
c
        return
        end
