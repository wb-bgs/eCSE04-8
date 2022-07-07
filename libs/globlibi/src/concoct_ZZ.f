cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine concoct_ZZ
c		V. Lesur 18/08/2011
c
c       Update a vector ZZ (nl1) by adding the diagonal elements of
c            A2'.W2'.R.W2.A2.ds
c       where A2 is the condition  matrix (nl2 lines and nl1 columns)
c             W2 is the weight matrix (diag nl2)
c             R  is the Misfit matrix (diag nl2)
c
c   Input:
c       nub     data type
c       DS      Vector of current descent direction (nl1)
c       HJ      Esimated diagonal of the Hessian (nl1)
c       A2	Matrice of equations of condition (nl2Xnl1)
c       W2s     Inverse of Diagonal covariance martix (W2 matrix squared) (nl2)
c       FM      Misfit Function 
c		   (user provided, external in calling program)
c       MV      Misfit to data vector (nl2)
c       nl1/2   number of lignes
c
c   output:
c	ZZ      UPDATED VECTOR 	AtWADS
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine concoct_ZZ(nub,FM,nl1,nl2,MV,W2s,A2,DS,ZZ)
c
        implicit none
c
        integer nl1,nl2,il,j,nub(*)
        real*8 DS(*),MV(*),W2s(*),A2(nl2,*),ZZ(*),dw1
c
        real*8 FM
        external FM
c
        do j=1,nl2
          dw1=dot_product(A2(j,1:nl1),DS(1:nl1))
          dw1=2.d0*FM(j,nub,MV)*W2s(j)*dw1
          do il=1,nl1
            ZZ(il)=ZZ(il)+dw1*A2(j,il)
          enddo
        enddo
c
        return
        end
