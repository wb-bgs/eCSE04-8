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
c       nub     data type
c       GJ      Weighted sum of square gradient (nl1)
c       A2      Matrice of equations of condition (nl2Xnl1)
c       W2s     Inverse of Diagonal covariance martix (W2 matrix squared) (nl2)
c       FM      Misfit Function 
c               (user provided, external in calling program)
c       MV      Misfit to data vector (nl2)
c       BB      [A2.M-B] (nl2)
c       nl1/2   number of lignes
c
c   output:
c       GJ      Updated gradient
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine concoct_GJ(nub,FM,nl1,nl2,MV,W2s,A2,BB,GJ)
c
        implicit none
c
        integer nl1,nl2,il,j,nub(*)
        real*8 MV(*),W2s(*),A2(nl2,*),BB(*),GJ(*),dw1
c
        real*8 FM
        external FM
c
        do j=1,nl2
          dw1=2.d0*FM(j,nub,MV)*W2s(j)*BB(j)
          do il=1,nl1
            GJ(il)=GJ(il)+A2(j,il)*dw1
          enddo
        enddo
c
        return
        end
