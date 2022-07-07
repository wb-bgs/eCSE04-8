C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
      subroutine gd2gc (cltgd,alt,xgd,zgd,cltgc,r,xgc,zgc)
c
c     conversion from geodetic to geocentric coordinates
c
c     IN: real*8
c      cltgd	geodetic colatitude (deg)
c      alt	altitude above sea-level (km)
c      xgd	north magnetic component in geodetic coords
c      zgd	vertical magnetic component in geodetic coords
c
c
c     OUT: real*8
c      cltgc	geocentic colatitude (deg)
c      r	radius from centre of earth (km)
c      xgc	north magnetic component in geocentric coords
c      zgc	vertical magnetic component in geocentric coords
c
c     Sep 2004 spheroid parameters changed from IAU66 to WGS84
c
      implicit double precision (a-h,o-z)
      one   = cltgd
      ct    = cosd(one)
      st    = sind(one)
c      a2    = 40680925.
c      b2    = 40408585.
      a2 = 40680631.6
      b2 = 40408296.0
      one   = a2*st*st
      two   = b2*ct*ct
      three = one + two
      rho   = sqrt(three)
      r     = sqrt(alt*(alt + 2.0*rho) + (a2*one + b2*two)/three)
      cd    = (alt + rho)/r
      sd    = (a2 - b2)/rho*ct*st/r
      one   = ct
      ct    = ct*cd -  st*sd
      st    = st*cd + one*sd
      cltgc = atan2d(st,ct)
      xgc   = xgd*cd - zgd*sd
      zgc   = zgd*cd + xgd*sd
c
      return
      end
