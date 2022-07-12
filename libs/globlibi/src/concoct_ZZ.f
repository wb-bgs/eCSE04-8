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
c       ds      Vector of current descent direction (nl1)
c       hj      Esimated diagonal of the Hessian (nl1)
c       a2	Matrice of equations of condition (nl2Xnl1)
c       w2s     Inverse of Diagonal covariance martix (W2 matrix squared) (nl2)
c       fun_mf  Misfit Function 
c		(user provided, external in calling program)
c       mv      Misfit to data vector (nl2)
c       nl1/2   number of lignes
c
c   output:
c	zz      UPDATED VECTOR 	AtWADS
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine concoct_ZZ(nub,fun_mf,nl1,nl2,mv,w2s,a2,ds,zz)
c
        implicit none
c
        integer nl1,nl2,il,j,nub(*)
        real*8 ds(*),mv(*),w2s(*),a2(nl2,*),zz(*),dw1
c
        real*8 fun_mf
        external fun_mf
c
        do j=1,nl2
          dw1=dot_product(a2(j,1:nl1),ds(1:nl1))
          dw1=2.d0*fun_mf(j,nub,mv)*w2s(j)*dw1
          do il=1,nl1
            zz(il)=zz(il)+dw1*a2(j,il)
          enddo
        enddo
c
        return
        end