cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine mklf_F
c		V. Lesur 01/02/2001
c
c    Modified 23/02/04 : GOTO removed
c
c       computes the Schmidt normalised legendre functions dlf for 
c       a given order "m" from degree m to lm (lm > m). 
c       Fast version: minimium check, do not avoid underflows
c       
c       recurrence 3.7.28 Fundations of geomagetism, Backus 1996
c
c    inputs:
c	dc/ds           cos(colatitude)/sin(colatitude)
c       nm/ilg          order and degree max 
c    outputs:
c       dlf             Legendre functions [dim min: max (2,lm-nm+1)]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine  mklf_F(nm, ilg, ds, dc, dlf)
c
        implicit none
c
        integer nm, ilg
        real*8 ds, dc, dlf(*)
c
        integer ierr, il, d0
        real*8 dnm, d1, d2, dbeta, dalpha
c
        real*8 dgamln
        external dgamln
c
c    start calculs
c
c       l=nm,  l=nm+1
c
        dnm = dble(nm)                    ! dble real for nm
        d1 = dgamln(2*dnm+1.0d0,ierr)     ! d1=log(fact(2dnm))
        d2 = dgamln(dnm+1.0d0,ierr)       ! d2=log(fact(dnm))
        if (ierr.ne.0) then
          write(*,*)'mklf_F: Cannot computes normalisation cst !'
          stop
        endif
c
        d2 = 0.5d0*d1 - d2                ! d2=sqrt(fact(2dnm))/fact(dnm)
        d2 = d2 - nm*dlog(2.0d0)          !
        d2 = dexp(d2)                     ! normalisation cst.
        if (nm.ne.0) d2 = d2*dsqrt(2.0d0) ! special case  m=0
c
        d1 = ds
        if (d1.ne.0.0d0) then
          d1 = d1**nm
        else
          if (nm.eq.0) d1 = 1.d0
        endif
c
        dlf(1) = d1*d2                              ! leg. func. (m,m)
        dlf(2) = dlf(1) * dc * dsqrt(dble(2*nm+1))  ! leg. func. (m+1,m)
c
c     l=nm+2....lm
c
        do il=2,ilg-nm                        ! l=nm+il-1
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