ccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Function l2_std
c                  V. Lesur	27/04/2005  
c
c     modified on 21.08.2012 such that std & npt are not modified 
c     during a call (V. Lesur)
c
c       Note: previously called l2_std_p. From 14.07.2011 replace the
c       previous version of l2_std (20.02.2002)
c
c       14.07.2011  nub added to parameter list
c
c      depending on input do one of the two:
c      1-  np.ge.0 
c            l2_std = std for np points added to a total of npt points
c      3-  np.lt.0 
c            l2_std = std of the combined (-np) stds     
c                   
c       All stds are computed with an L2 norm
c
c     input np>0:
c       nub     type of data
c       npt     np total number of point used to compute std
c       std     values of std computed with npt points
c       np      number of point added in this itteration
c	mv	misfit vector dimmin np
c       wgh	vector of 1/variances dimmin np
c     input np<0:
c       nub     type of data
c       np      (-1*) number of std to be added
c	mv	std values  dimmin -np
c       wgh	number of point used to compute stds dimmin -np
c    
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 function l2_std(nub,npt,np,std,mv,wgh) 
c
        implicit none
c
        integer npt,np,i,nub(*),nn
        real*8 mv(*),std,wgh(*),dd
c
        if(np.ge.0)then
          dd=dble(npt)*std**2
          do i=1,np
            dd=dd+wgh(i)*mv(i)**2
          enddo
c
          l2_std=dsqrt(dd/dble(npt+np))
        else
          dd=0.0d0
          nn=0
          do i=1,-np
            dd=dd+wgh(i)*mv(i)**2
            nn=nn+idnint(wgh(i))
          enddo
c
          l2_std=dsqrt(dd/dble(nn))
        endif
c
        return
        end
