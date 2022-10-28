cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine mklf_F2
c		M. Bareford 28/10/2022
c
c       Identical to  mklf_F() but with a precomputed d2a array:
c       this avoids having to repeatedly call dgamln.     
c
c    inputs:
c       nm/ilg          order and degree max
c	dc/ds           cos(colatitude)/sin(colatitude)
c       d2a             pre-computed array 
c    outputs:
c       dlf             Legendre functions [dim min: max (2,lm-nm+1)]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine  mklf_F2(nm, ilg, ds, dc, d2a, dlf)
c
        implicit none
c
        integer nm, ilg
        real*8 ds, dc, d2a(0:ilg), dlf(*)
c
        integer il, d0
        real*8 d1, d2, dbeta, dalpha
c
        d1 = ds
        if (d1.ne.0.0d0) then
          d1 = d1**nm
        else
          if (nm.eq.0) d1 = 1.d0
        endif
c
        dlf(1) = d1*d2a(nm)                         ! leg. func. (m,m)
        dlf(2) = dlf(1) * dc * dsqrt(dble(2*nm+1))  ! leg. func. (m+1,m)
c
c     l=nm+2....lm
c
        do il=2,ilg-nm                       ! l=nm+il-1
c
          d0 = il+2*nm
          d1 = dble((il-1) * (d0-1))         ! (l-m)*(l+m)
          d2 = dble(il * d0)                 ! (l-m+1)*(l+m+1)	
          dbeta = dsqrt(d1/d2)               ! recurrence coeff.
          d1 = dble(2*(il+nm)-1)             ! 2l+1
          dalpha = d1/dsqrt(d2)              !
c
          ! leg. func. (nm+il-1,nm) 
          dlf(il+1) = dalpha*dlf(il)*dc - dbeta*dlf(il-1)

        enddo
c
        return
        end