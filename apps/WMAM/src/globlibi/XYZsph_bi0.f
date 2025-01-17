        module XYZsph_bi0
c
          implicit none
c
          public
     >      XYZsph_bi0_sample,
     >      XYZsph_bi0_fun,
     >      XYZsph_bi0_sub

          private
     >      mk_lf_dlf
c
c
        contains
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine XYZsph_bi0_sample
c		
c       Compute dx = bx.bc, dy = by.bc and dz = bz.bc.
c       This subroutine is called for sample points only.
c             
c       input:
c          shdeg  INTEGER       max SH degree value
c          nb     INTEGER       number of coefficients
c          d2a    REAL*8        pre-computed array used by mk_lf_dlf() subroutine
c          (d)dlf REAL*8        pre-allocated arrays computed by mk_lf_dlf() and
c                               used within this subroutine
c          bc     REAL*8        coefficient array
c          p1     REAL*8        co-latitude
c          p2     REAL*8        longitude
c          ra     REAL*8        radius
c
c       output:
c          bex    REAL*8        x component of magnetic field
c          bey    REAL*8        y component of magnetic field
c          bez    REAL*8        z component of magnetic field
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        attributes(device)
     >  subroutine XYZsph_bi0_sample(shdeg, nb, d2a,
     >                               dlf, ddlf, bc,
     >                               p1, p2, ra,
     >                               bex, bey, bez)
c
        implicit none
c
        integer shdeg, nb
        real(8) d2a(0:shdeg)
        real(8) dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real(8) bc(1:nb)
        real(8) p1, p2, ra
        real(8) bex, bey, bez 
c
        integer nu, il, im, ik
        real(8) rc, rs
        real(8) ds, dc, dr, dw
        real(8) bx, by, bz
        real(8) bxp1, byp1, bzp1
        real(8) bc_nu, bc_nup1
c
        real(8) dx, dy, dz, dd
        real(8) dxbey, dxbez
        real(8) dybex, dybez
        real(8) dzbex, dzbey
c
        real(8) xy_c, xz_c
        real(8) yx_c, yz_c
        real(8) zx_c, zy_c
c
        real(8) bex2, bey2, bez2
c        
c 
        dx = 0.0d0
        dy = 0.0d0 
        dz = 0.0d0 
c
        rc = dcos(p1)
        rs = dsin(p1)
c
        call mk_lf_dlf(0, shdeg, rs, rc,
     >                 d2a(0), dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
c          
          ik = il+1
c          
          bx =  ddlf(ik) * dr
          bz = -dlf(ik)  * dr * dble(ik)
c
          dx = dx + bx * bc(il)
          dz = dz + bz * bc(il)
        enddo
c
c
        nu = shdeg + 1
        do im = 1,shdeg
          call mk_lf_dlf(im, shdeg, rs, rc,
     >                   d2a(im), dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            bc_nu   = bc(nu)
            bc_nup1 = bc(nu+1)

            ik = il-im+1
            dr = ra**(il+2)

            dw   = ddlf(ik) * dr
            bx   = dc * bc_nu
            bxp1 = ds * bc_nup1
            dx   = dx + dw*(bx + bxp1)

            dw   = (dlf(ik)/rs) * dr * dble(im)
            by   =  ds * bc_nu
            byp1 = -dc * bc_nup1
            dy   =  dy + dw*(by + byp1)
            
            dw   =  dlf(ik) * dr * dble(il+1)
            bz   = -dc * bc_nu
            bzp1 = -ds * bc_nup1
            dz   =  dz + dw*(bz + bzp1)

            nu = nu + 2
          enddo
        enddo
c
        dxbey = dx*bey
        dxbez = dx*bez
        dybex = dy*bex
        dybez = dy*bez
        dzbex = dz*bex
        dzbey = dz*bey
c
        xy_c = dxbey - dybex
        xz_c = dxbez - dzbex
        yx_c = -xy_c
        yz_c = dybez - dzbey
        zx_c = -xz_c
        zy_c = -yz_c
c
        dd = dsqrt(yz_c**2 + xz_c**2 + xy_c**2)
c
        bex2 = (xz_c*bez + xy_c*bey) / dd
        bey2 = (yz_c*bez + yx_c*bex) / dd
        bez2 = (zy_c*bey + zx_c*bex) / dd
c	
        bex = bex2
        bey = bey2
        bez = bez2
c
        end subroutine XYZsph_bi0_sample
c
c     
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	function XYZsph_bi0_fun
c		
c       Computes the dot product of 'be' and 'bc', avoiding the need
c       to store the entire contents of the 'be' coefficient array.
c             
c       input:
c          shdeg  INTEGER     max SH degree value
c          nb     INTEGER     number of coefficients
c          d2a    REAL*8      pre-computed array used by mk_lf_dlf() subroutine
c          (d)dlf REAL*8      pre-allocated arrays computed by mk_lf_dlf() and
c                             used within this function
c          bc     REAL*8      coefficient array
c          p1     REAL*8      co-latitude
c          p2     REAL*8      longitude
c          ra     REAL*8      radius
c          bex    REAL*8      x component of magnetic field
c          bey    REAL*8      y component of magnetic field
c          bez    REAL*8      z component of magnetic field
c
c       output:
c          YZsph_bi0_fun  REAL*8
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        attributes(device)
     >  real*8 function XYZsph_bi0_fun(shdeg, nb, d2a,
     >                                 dlf, ddlf, bc,
     >                                 p1, p2, ra,
     >                                 bex, bey, bez)
c
        implicit none
c
        integer shdeg, nb
        real(8) d2a(0:shdeg)
        real(8) dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real(8) bc(1:nb)
        real(8) p1, p2, ra
        real(8) bex, bey, bez 
c
        integer nu, il, im, ik
        real(8) rc, rs
        real(8) ds, dc, dr, dw
        real(8) bx, by, bz
        real(8) bxp1, byp1, bzp1
c 
c
        XYZsph_bi0_fun = 0.0d0 
c        
        rc = dcos(p1)
        rs = dsin(p1)
c
        call mk_lf_dlf(0, shdeg, rs, rc,
     >                 d2a(0), dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
c          
          ik = il+1
c          
          bx =  ddlf(ik) * dr
          bz = -dlf(ik)  * dr * dble(ik)
c
          XYZsph_bi0_fun = XYZsph_bi0_fun
     >                   + (bex*bx + bez*bz) * bc(il)
        enddo
c
c
        nu = shdeg + 1
        do im = 1,shdeg
          call mk_lf_dlf(im, shdeg, rs, rc,
     >                   d2a(im), dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            ik = il-im+1
            dr = ra**(il+2)
c
            dw   = ddlf(ik) * dr * bex
            bx   = dw*dc
            bxp1 = dw*ds
c            
            dw   = (dlf(ik)/rs) * dr * dble(im) * bey
            by   =  dw*ds
            byp1 = -dw*dc
c            
            dw   =  dlf(ik) * dr * dble(il+1) * bez
            bz   = -dw*dc
            bzp1 = -dw*ds
c
            XYZsph_bi0_fun = XYZsph_bi0_fun
     >                     + (bx + by + bz) * bc(nu)
     >                     + (bxp1 + byp1 + bzp1) * bc(nu+1)

            nu = nu + 2
          enddo
        enddo
c
        end function XYZsph_bi0_fun
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine XYZsph_bi0_sub
c		
c       Computes the gradient of the weighted sum of squares (gj) and
c       the diagonal of the Hessian (hj). 
c
c       input:
c          shdeg    INTEGER     max SH degree value
c          nb       INTEGER     number of coefficients
c          d2a      REAL*8      pre-computed array used by mk_lf_dlf() subroutine
c          (d)dlf   REAL*8      pre-allocated arrays computed by mk_lf_dlf() and
c                               used within this subroutine
c          p1       REAL*8      co-latitude
c          p2       REAL*8      longitude
c          ra       REAL*8      radius
c          bex      REAL*8      x component of magnetic field
c          bey      REAL*8      y component of magnetic field
c          bez      REAL*8      z component of magnetic field 
c          dw_gj    REAL*8      multiplier for gj terms
c          dw_dh    REAL*8      multiplier for dh terms
c
c       output:
c          gj       REAL*8      gradient of the weighted sum of squares (nb)
c          dh       REAL*8      diagonal of the Hessian (nb)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        attributes(device)
     >  subroutine XYZsph_bi0_sub(shdeg, nb, d2a,
     >                            dlf, ddlf,
     >                            p1, p2, ra,
     >                            bex, bey, bez,
     >                            dw_gj, dw_dh,
     >                            gj, dh)
c
        use cudafor
c
        implicit none
c
        integer shdeg, nb
        real(8) d2a(0:shdeg)
        real(8) dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real(8) p1, p2, ra
        real(8) bex, bey, bez 
        real(8) dw_gj, dw_dh
        real(8) gj(1:nb), dh(1:nb) 
c
        integer nu, il, im, ik
        real(8) rc, rs
        real(8) ds, dc, dr, dw
        real(8) bx, by, bz
        real(8) bxp1, byp1, bzp1
        real(8) be, bep1
        real(8) oldval
c 
c
        rc = dcos(p1)
        rs = dsin(p1)
c
        call mk_lf_dlf(0, shdeg, rs, rc,
     >                 d2a(0), dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
c          
          ik = il+1
c          
          bx =  ddlf(ik) * dr
          bz = -dlf(ik)  * dr * dble(ik)
c
          be = (bex*bx + bez*bz)
c
          oldval = atomicadd(gj(il), dw_gj*be)
          oldval = atomicadd(dh(il), dw_dh*(be**2))
c
        enddo
c
c
        nu = shdeg + 1
        do im = 1,shdeg
          call mk_lf_dlf(im, shdeg, rs, rc,
     >                   d2a(im), dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            ik = il-im+1
            dr = ra**(il+2)
c
            dw   = ddlf(ik) * dr
            bx   = dw*dc
            bxp1 = dw*ds
c            
            dw   = (dlf(ik)/rs) * dr * dble(im)
            by   =  dw*ds
            byp1 = -dw*dc
c            
            dw   =  dlf(ik) * dr * dble(il+1)
            bz   = -dw*dc
            bzp1 = -dw*ds
c
            be   = bex*bx + bey*by + bez*bz
            bep1 = bex*bxp1 + bey*byp1 + bez*bzp1
c
            oldval = atomicadd(gj(nu), dw_gj*be)
            oldval = atomicadd(gj(nu+1), dw_gj*bep1)
c
            oldval = atomicadd(dh(nu), dw_dh*(be**2))
            oldval = atomicadd(dh(nu+1), dw_dh*(bep1**2))
c
            nu = nu + 2
          enddo
        enddo
c
        end subroutine XYZsph_bi0_sub
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 	subroutine mk_lf_dlf
c               V. Lesur 19 June 2005
c
c       Computes legendre Functions and their derivatives along theta
c
c       CALLED: mklf_F.f
c
c       limitations of mklf_F apply
c
c       Fast version: tests reduced to minimum
c       test in loop in else branch eliminated  
c
c       input:
c         im/shdeg      order and degree max
c         rc/rs         cos(colatitude)/sin(colatitude)
c         d2a_im        element at position im from d2a array
c       output:
c         dlf           legendre function from im to nl
c         ddlf          derivative of legendre function from im to nl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        attributes(device)
     >  subroutine mk_lf_dlf(im, shdeg, rs, rc, 
     >                       d2a_im, dlf, ddlf)
c
        implicit none
c
        integer im, shdeg
        real(8) rs, rc, d2a_im
        real(8) dlf(1:shdeg+1), ddlf(1:shdeg+1)
c
        integer d0, il, jl
        real(8) d1, d2
        real(8) dalpha, dbeta
c
c
c  Initialise dlf array
        d1 = rs
        if (d1 .ne. 0.0d0) then
          d1 = d1**im
        else
          if (im .eq. 0) d1 = 1.d0
        endif
c
        dlf(1) = d1*d2a_im
        dlf(2) = dlf(1) * rc * dsqrt(dble(2*im+1))
c
        do il = 2,shdeg-im
          d0 = il+2*im
          d1 = dble((il-1) * (d0-1))
          d2 = dble(il * d0)
          dbeta = dsqrt(d1/d2)*dlf(il-1)
c
          d1 = dble(2*(il+im)-1)
          dalpha = (d1/dsqrt(d2))*dlf(il)*rc
c
          dlf(il+1) = dalpha - dbeta
        enddo
c
c
c  Initialise ddlf array
        if (rs .eq. 0.0d0) then
c
          if (im .ne. 1) then 
            ddlf(1:shdeg+1) = 0.0d0
          else
            ddlf(1) = -1.0d0
            ddlf(2) = -dsqrt(3.0d0) 
            do il = 3,shdeg+1
              jl = (il**2)-1
              d1 = (2*il-1)/dsqrt(dble(jl))
              d2 = dsqrt(dble((il-1)**2-1)/jl)
              ddlf(il) = d1*ddlf(il-1) - d2*ddlf(il-2)
            enddo
          endif
c
        else
c
          do il = shdeg,im+1,-1
            jl = il - im + 1
            d1 = dsqrt(dble((il-im)*(il+im)))
            d2 = dble(il) 
            ddlf(jl) = (d2*rc*dlf(jl) - d1*dlf(jl-1)) / rs
          enddo
          ddlf(1) = im * rc * dlf(1)/rs
c
        endif
c
        end
c
c
        end module XYZsph_bi0