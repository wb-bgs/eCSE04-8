cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine mklf_F2
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
c         ilg  INTEGER  degree max
c         nm   INTEGER  order
c         ds   REAL*8  sin(colatitude)
c         dc   REAL*8  cos(colatitude)
c         d2a, d3a, dalphaa, dbetaa  REAL*8
c           arrays holding pre-computed values for mklf_F2()
c    outputs:
c       dlf  REAL*8  Legendre functions [dim min: max (2,lm-nm+1)]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine  mklf_F2(nm, ilg, d2a, d3a,
     >                      dalphaa, dbetaa,
     >                      ds, dc, dlf)
c
        implicit none
c
        integer nm, ilg
        real*8 d2a(0:ilg), d3a(0:ilg)
        real*8 dalphaa(*), dbetaa(*)
        real*8 ds, dc
        real*8 dlf(*)
c
        integer il
        real*8 d1
        integer, save :: k=0
c
c
        d1 = ds
        if (d1.ne.0.0d0) then
          d1 = d1**nm
        else
          if (nm.eq.0) d1 = 1.d0
        endif
c
        dlf(1) = d1*d2a(nm)
        dlf(2) = dlf(1) * dc * d3a(nm)
c
        if (nm .eq. 0) then
          k = 1
        endif
c
        do il = 2,ilg-nm
          dlf(il+1) = dalphaa(k)*dlf(il)*dc - dbetaa(k)*dlf(il-1)
          k = k+1
        enddo
c
        return
        end