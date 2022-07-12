cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine concoct_HJ
c		V. Lesur 08/06/2006
c
c       14.07.2011  nub added to parameter list
c                   implicit none
c
c       Update a vector HJ (nl1) by adding the diagonal elements of
c            A2'.W2'.R.W2.A2
c       where A2 is the condition  matrix (nl2 lines and nl1 columns)
c             W2 is the weight matrix (diag nl2)
c             R  is the Misfit matrix (diag nl2)
c
c   Input:
c       nub     data type
c       hj      Esimated diagonal of the Hessian (nl1)
c       a2	Matrice of equations of condition (nl2Xnl1)
c       w2s     Inverse of Diagonal covariance martix (W2 matrix squared) (nl2)
c       fun_mf  Misfit Function 
c		(user provided, external in calling program)
c       mv      Misfit to data vector (nl2)
c       nl1/2   number of lignes
c
c   output:
c	hj      Updated diagonal of Hessian
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine concoct_HJ(nub,fun_mf,nl1,nl2,mv,w2s,a2,hj)
c
        implicit none
c
        integer nl1,nl2,il,j,nub(*)
        real*8 mv(*),w2s(*),a2(nl2,*),hj(*),dw1
c
        real*8 fun_mf
        external fun_mf
c
        do j=1,nl2
          dw1=2.d0*fun_mf(j,nub,mv)*w2s(j)
          do il=1,nl1
            hj(il)=hj(il)+dw1*a2(j,il)**2
          enddo
        enddo
c
        return
        end