c====================================================================
c
c   *** GAULEG - find zeroes and weights for Gauss Quadrature ***
c     probably a Gubbins subroutine originally
c     check 31/05/2004 up to n=96 against  Abrm&Stegun p. 919: 
c                      GOTO removed
c                      error on starting value z corrected 
c                      maximum number of iteration introduced
c   
c           
c
c   PARAMETERS
c        x1,x2: end points of integration
c        x/w:   arrays of zeroes of Legendre function of degree n
c               corresponding weights.
c        n:     degree of integration
c
c====================================================================
        subroutine gauleg(x1,x2,x,w,n)
c
        implicit real*4 (a-c,e-h,o-z)
        implicit real*8 (d)
        implicit integer*4 (i-n)
c
        parameter (itmax=10)
c
        real*8 x(n), w(n), dpi, deps
        real*8 x1, x2, xm, xl
        real*8 p1, p2, p3, pp, z, z1
c
        deps = 1.d-15
        dpi = 4.d0*datan(1.d0)
c
        m = int(float(n+1)/2.)
        x(m+1) = 0.0d0
c
        xm = 0.5d0*(x2+x1)
        xl = 0.5d0*(x2-x1)
c
        do i = 1,m
c
c  Search x value by gradient methods
c
            z = dcos(dpi*dble(2*i)/dble(2*n+1))
            z1 = dcos(dpi*dble(2*i-1)/dble(2*n+1))
            it = 0
            do while (dabs(z-z1) .gt. deps 
     >                .and.
     >                it .le. itmax)
                it = it+1
c
c  Recurence on Pl(cos(z))
                p1 = 1.d0
                p2 = 0.0d0
                do j = 1,n
                    p3 = p2
                    p2 = p1
                    p1 = (2.d0*j-1.d0)*z*p2
     >                   - (dble(j)-1.d0)*p3
                    p1 = p1/dble(j)
                enddo
c
c  Calculate derivative relative to cos(z)
                pp = n*(z*p1-p2)/(z*z-1.d0)
c
c  Newton's iteration
                z1 = z
                z = z1-p1/pp
            enddo
c
           if (it .ge. itmax) then
               write(*,*)'Warning: convergence for zero:',i
           endif
c
c  Found saved + weight calculated (twice by symmetry)
           x(i) = xm-xl*z
           x(n+1-i) = xm+xl*z
           w(i) = xl*2.d0/((1.d0-z*z)*pp*pp)
           w(n+1-i) = w(i)
c
        enddo
c
        return
        end