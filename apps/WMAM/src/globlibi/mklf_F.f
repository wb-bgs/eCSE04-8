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
        real*8 ds, dc, dlf(1:ilg+1)
c
        integer ierr, il, d0
        real*8 dnm, d1, d2, dbeta, dalpha
c
        real*8 dgamln
        external dgamln
c
c
        dnm = dble(nm)
        d1 = dgamln(2*dnm+1.0d0,ierr)
        d2 = dgamln(dnm+1.0d0,ierr)
        if (ierr .ne. 0) then
          write(*,*) 'mklf_F: Cannot computes normalisation cst !'
          stop
        endif
c
        d2 = 0.5d0*d1 - d2
        d2 = d2 - nm*dlog(2.0d0)
        d2 = dexp(d2)
        if (nm .ne. 0) then
          d2 = d2*dsqrt(2.0d0)
        endif
c
        d1 = ds
        if (d1 .ne. 0.0d0) then
          d1 = d1**nm
        else
          if (nm .eq. 0) d1 = 1.d0
        endif
c
        dlf(1) = d1*d2
        dlf(2) = dlf(1) * dc * dsqrt(dble(2*nm+1))
c
        do il = 2,ilg-nm
c
          d0 = il + 2*nm
          d1 = dble((il-1) * (d0-1))
          d2 = dble(il * d0)
          dbeta = dsqrt(d1/d2)
          d1 = dble(2*(il+nm)-1)
          dalpha = d1/dsqrt(d2)
c
          dlf(il+1) = dalpha*dlf(il)*dc - dbeta*dlf(il-1)
c
        enddo
c
        return
        end