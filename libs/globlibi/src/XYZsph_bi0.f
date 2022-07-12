cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Subroutine XYZsph_bi0
c		Vincent Lesur 21/01/2009
c
c    modified:
c       07.01.2010  V.Lesur :
c         1- reference radius ra as parameter
c
c       For a all spherical harmonique between degree 1 and nlie
c       and a  colatitude, longitude, radius computes for
c       internal sources: 
c       X component basis  (X horizontal North)
c       Y component basis  (Y horizontal East)
c       Z component basis  (Z Vertical Down) 
c
c       These basis is the derivative of a potential 
c       defined on the sphere of radius "ra"  
c       (See Foundations of Geomagnetism, Backus 1996 page 110 and 125)
c      
c       No time dependence
c
c       Spherical harmonics are set by increasing degree "l", for 
c       a given degree by increasing order "m" (l=0,m=0 excluded)
c       nu=(l**2+2*(m-1))   FOR cosine terms if m.ne.0
c       nu=(l**2+2*m-1)     For sine terms or m=0
c
c       colat,long,radius given in degree,degree,km
c
c       input:
c         ilg  INTEGER   max SH degree 
c         rag  REAL*8    reference radius
c         d2a, d3a, dalphaa, dbetaa  REAL*8
c           arrays holding pre-computed values for mklf_F2()
c         d4a, d5a, d6a  REAL*8
c           arrays holding pre-computed values for mk_lf_dlf()
c         pos  REAL*8    positions in space and time (colat,long,radius)
c
c       output:
c         bx   REAL*8   X Base funtion value
c         by   REAL*8   Y Base funtion value
c         bz   REAL*8   Z Base funtion value
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine XYZsph_bi0(ilg, rag, d2a, d3a,
     >                        dalphaa, dbetaa,
     >                        d4a, d5a, d6a,
     >                        pos, bx, by, bz)
c
        implicit none
c
        integer ilg
        real*8 rag
        real*8 d2a(0:ilg), d3a(0:ilg)
        real*8 dalphaa(*), dbetaa(*)
        real*8 d4a(1:ilg+1), d5a(3:ilg+1), d6a(*)
        real*8 pos(*), bx(*), by(*), bz(*)
c
        integer il, im, nu
        real*8 rc, rs, dw
        real*8 ds, dc, dcosd, dsind, ra_div_pos3
        real*8, allocatable :: dlf(:), ddlf(:), dr(:)
c
        external dcosd,dsind
c
        allocate(dlf(1:ilg+1), ddlf(1:ilg+1))
        allocate(dr(1:ilg))
c
        rc = dcosd(pos(1))
        rs = dsind(pos(1))

        ra_div_pos3 = rag / pos(3)
c
c   im=0
        im=0
        call mk_lf_dlf(im, ilg, d2a, d3a,
     >                 dalphaa, dbetaa,
     >                 d4a, d5a, d6a,
     >                 rs, rc, dlf, ddlf)

        do il=1,ilg
          nu = il*il - 1
          dr(il) = ra_div_pos3**(il+2)
          bx(nu+1) = ddlf(il-im+1) * dr(il)
          by(nu+1) = 0.0d0
          bz(nu+1) = -dlf(il-im+1) * dble(il+1) * dr(il)
        enddo

c   im.ne.0
        do im=1,ilg
          call mk_lf_dlf(im, ilg, d2a, d3a,
     >                   dalphaa, dbetaa,
     >                   d4a, d5a, d6a,
     >                   rs, rc, dlf, ddlf)

          dc = dcosd(im*pos(2))
          ds = dsind(im*pos(2))
          do il=im,ilg
            nu = il*il + 2*(im-1)

            dw = ddlf(il-im+1) * dr(il)
            bx(nu+1) = dw*dc
            bx(nu+2) = dw*ds
            
            dw = dlf(il-im+1)/rs
            dw = dw * dble(im) * dr(il)
            by(nu+1) = dw*ds
            by(nu+2) = -dw*dc
            
            dw = dlf(il-im+1) * dble(il+1) * dr(il)
            bz(nu+1) = -dw*dc
            bz(nu+2) = -dw*ds
          enddo
        enddo

        deallocate(dlf, ddlf, dr)
c
        return
        end