cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine concoct_GJ
c		V. Lesur 11/02/2005
c
c       14.07.2011  nub added to parameter list
c
c       Update a vector GJ (nl1) by adding the product
c            A2'.W2'.R.W2.[A2.M-B]
c       where A2 is the condition  matrix (nl2 lines and nl1 columns)
c             W2 is the weight matrix (diag nl2)
c             M  is the estimated parameter vector(nl1)
c             R  is the Misfit matrix (diag nl2)
c             B  data vector (nl2)
c
c   Input:
c       fun_mf  Misfit Function 
c               (user provided, external in calling program)
c       nl1/2   number of lignes
c       w2s     Inverse of Diagonal covariance martix (W2 matrix squared) (nl2)
c       a2      Matrice of equations of condition (nl2Xnl1)
c       bb      [A2.M-B] (nl2)
c       gj      Weighted sum of square gradient (nl1)
c
c   output:
c       gj      Updated gradient
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine concoct_GJ(fun_mf,nl1,nl2,w2s,a2,bb,gj)
c
        implicit none
c
        integer nl1,nl2,il,j
        real*8 w2s(*),a2(nl2,*),bb(*),gj(*),dw1
c
        real*8 fun_mf
        external fun_mf
c
        do j=1,nl2
          dw1=2.d0*fun_mf()*w2s(j)*bb(j)
          do il=1,nl1
            gj(il)=gj(il)+a2(j,il)*dw1
          enddo
        enddo
c
        return
        end