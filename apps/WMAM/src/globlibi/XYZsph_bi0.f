ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine XYZsph_bi0_sample
c		
c       Compute dx = bx.bc, dy = by.bc and dz = bz.bc.
c       This subroutine is called for sample points only.
c             
c       input:
c          shdeg  INTEGER       max SH degree value
c          nb     INTEGER       number of coefficients
c          d2a    REAL*8        pre-computed array for mk_lf_dlf() subroutine
c          dra    REAL*8        pre-allocated array used in this subroutine
c          dalpha REAL*8        pre-allocated array used in mk_lf_dlf() subroutine
c          dbeta  REAL*8        "
c          dlf    REAL*8        "
c          ddlf   REAL*8        "
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
        subroutine XYZsph_bi0_sample(shdeg, nb, d2a,
     >                               dra, dalpha, dbeta,
     >                               dlf, ddlf, bc,
     >                               p1, p2, ra,
     >                               bex, bey, bez)
c
        implicit none
c
        integer shdeg, nb
        real*8 d2a(0:shdeg), dra(1:shdeg)
        real*8 dalpha(0:shdeg), dbeta(0:shdeg)
        real*8 dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real*8 bc(1:nb)
        real*8 p1, p2, ra
        real*8 bex, bey, bez 
c
        integer nu, il, im, ik
        real*8 rc, rs
        real*8 ds, dc, dr, dw
        real*8 bx, by, bz
        real*8 bxp1, byp1, bzp1
        real*8 bc_nu, bc_nup1
c
        real*8 dx, dy, dz, dd
        real*8 dxbey, dxbez
        real*8 dybex, dybez
        real*8 dzbex, dzbey
c
        real*8 xy_c, xz_c
        real*8 yx_c, yz_c
        real*8 zx_c, zy_c
c
        real*8 bex2, bey2, bez2
c        
c 
#if defined(OMP_OFFLOAD_CPTP) || defined(OMP_OFFLOAD_SSQGH)
!$omp declare target
#endif
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
     >                 d2a(0), dalpha, dbeta,
     >                 dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
          dra(il) = dr
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
     >                   d2a(im), dalpha, dbeta,
     >                   dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            bc_nu   = bc(nu)
            bc_nup1 = bc(nu+1)

            ik = il-im+1
            dr = dra(il)

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
        return
        end subroutine XYZsph_bi0_sample

      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	function XYZsph_bi0_fun
c		
c       Computes the dot product of 'be' and 'bc', avoiding the need
c       to store the entire contents of the 'be' coefficient array.
c             
c       input:
c          shdeg  INTEGER     max SH degree value
c          nb     INTEGER     number of coefficients
c          d2a    REAL*8      pre-computed array for mk_lf_dlf() subroutine
c          dra    REAL*8      pre-allocated array used in this subroutine
c          dalpha REAL*8      pre-allocated array used in mk_lf_dlf() subroutine
c          dbeta  REAL*8      "
c          dlf    REAL*8      "
c          ddlf   REAL*8      "
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
        real*8 function XYZsph_bi0_fun(shdeg, nb, d2a,
     >                                 dra, dalpha, dbeta,
     >                                 dlf, ddlf, bc,
     >                                 p1, p2, ra,
     >                                 bex, bey, bez)
c
        implicit none
c
        integer shdeg, nb
        real*8 d2a(0:shdeg), dra(1:shdeg)
        real*8 dalpha(0:shdeg), dbeta(0:shdeg)
        real*8 dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real*8 bc(1:nb)
        real*8 p1, p2, ra
        real*8 bex, bey, bez 
c
        integer nu, il, im, ik
        real*8 rc, rs
        real*8 ds, dc, dr, dw
        real*8 bx, by, bz
        real*8 bxp1, byp1, bzp1
c 
c
#if defined(OMP_OFFLOAD_CPTP)
!$omp declare target
#endif
c
c
        XYZsph_bi0_fun = 0.0d0 
c        
        rc = dcos(p1)
        rs = dsin(p1)
c
        call mk_lf_dlf(0, shdeg, rs, rc,
     >                 d2a(0), dalpha, dbeta,
     >                 dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
          dra(il) = dr
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
     >                   d2a(im), dalpha, dbeta,
     >                   dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            ik = il-im+1
            dr = dra(il)
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
        return
        end function XYZsph_bi0_fun


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine XYZsph_bi0_sub
c		
c       Computes the gradient of the weighted sum of squares (gj) and
c       the diagonal of the Hessian (hj). 
c
c       input:
c          shdeg    INTEGER     max SH degree value
c          nb       INTEGER     number of coefficients
c          d2a      REAL*8      pre-computed array for mk_lf_dlf() subroutine
c          dra      REAL*8      pre-allocated array used in this subroutine
c          dalpha   REAL*8      pre-allocated array used in mk_lf_dlf() subroutine
c          dbeta    REAL*8      "
c          dlf      REAL*8      "
c          ddlf     REAL*8      "
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
        subroutine XYZsph_bi0_sub(shdeg, nb, d2a,
     >                            dra, dalpha, dbeta,
     >                            dlf, ddlf,
     >                            p1, p2, ra,
     >                            bex, bey, bez,
     >                            dw_gj, dw_dh,
     >                            gj, dh)
c
        implicit none
c
        integer shdeg, nb
        real*8 d2a(0:shdeg), dra(1:shdeg)
        real*8 dalpha(0:shdeg), dbeta(0:shdeg)
        real*8 dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real*8 p1, p2, ra
        real*8 bex, bey, bez 
        real*8 dw_gj, dw_dh
        real*8 gj(1:nb), dh(1:nb) 
c
        integer nu, il, im, ik
        real*8 rc, rs
        real*8 ds, dc, dr, dw
        real*8 bx, by, bz
        real*8 bxp1, byp1, bzp1
        real*8 be, bep1
c 
c
#if defined(OMP_OFFLOAD_SSQGH)
!$omp declare target
#endif
c
c
        rc = dcos(p1)
        rs = dsin(p1)
c
        call mk_lf_dlf(0, shdeg, rs, rc,
     >                 d2a(0), dalpha, dbeta,
     >                 dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
          dra(il) = dr
c          
          ik = il+1
c          
          bx =  ddlf(ik) * dr
          bz = -dlf(ik)  * dr * dble(ik)
c
          be = (bex*bx + bez*bz)
c
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP ATOMIC UPDATE
          gj(il) = gj(il) + dw_gj*be
!$OMP ATOMIC UPDATE
          dh(il) = dh(il) + dw_dh*(be**2)
#else
          gj(il) = gj(il) + dw_gj*be
          dh(il) = dh(il) + dw_dh*(be**2)
#endif
c
        enddo
c
c
        nu = shdeg + 1
        do im = 1,shdeg
          call mk_lf_dlf(im, shdeg, rs, rc,
     >                   d2a(im), dalpha, dbeta,
     >                   dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            ik = il-im+1
            dr = dra(il)
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
#if defined(OMP_OFFLOAD_SSQGH)
!$OMP ATOMIC UPDATE
            gj(nu)   = gj(nu)   + dw_gj*be
!$OMP ATOMIC UPDATE
            gj(nu+1) = gj(nu+1) + dw_gj*bep1
!$OMP ATOMIC UPDATE
            dh(nu)   = dh(nu)   + dw_dh*(be**2)
!$OMP ATOMIC UPDATE
            dh(nu+1) = dh(nu+1) + dw_dh*(bep1**2)
#else
            gj(nu)   = gj(nu)   + dw_gj*be
            gj(nu+1) = gj(nu+1) + dw_gj*bep1
            dh(nu)   = dh(nu)   + dw_dh*(be**2)
            dh(nu+1) = dh(nu+1) + dw_dh*(bep1**2)
#endif
c
            nu = nu + 2
          enddo
        enddo
c
        return
        end subroutine XYZsph_bi0_sub