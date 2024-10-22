cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine XYZsph_bi0_sample
c		
c       Compute dx = bx.bc, dy = by.bc and dz = bz.bc.
c       This subroutine is called for sample points only.
c             
c       input:
c          shdeg  INTEGER       max SH degree value
c          nb     INTEGER       number of coefficients
c          d2a    REAL*8        pre-computed array for mk_lf_dlf()
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
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine XYZsph_bi0_sample(shdeg, nb,
     >                               d2a, bc,
     >                               p1, p2, ra,
     >                               bex, bey, bez)
c
        implicit none
c
        integer shdeg, nb
        real*8 d2a(0:shdeg), bc(nb)
        real*8 p1, p2, ra
        real*8 bex, bey, bez 
c
        integer nu, il, im, ik
        real*8 rc, rs
        real*8 ds, dc, dr, dw
        real*8 bx, by, bz
        real*8 bxp1, byp1, bzp1
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
        real*8, allocatable :: dalpha(:), dbeta(:)
        real*8, allocatable :: dra(:), dlf(:), ddlf(:)
c        
c 
#if defined(OMP_OFFLOAD_CPTP) || defined(OMP_OFFLOAD_SSQGH)
!$omp declare target
#endif
c
c
        allocate(dalpha(2:shdeg-1), dbeta(2:shdeg-1))
        allocate(dra(shdeg), dlf(shdeg+1), ddlf(shdeg+1))
c
        dx = 0.0d0
        dy = 0.0d0 
        dz = 0.0d0 
c
        rc = dcos(p1)
        rs = dsin(p1)
c
        im=0
        call mk_lf_dlf(im, shdeg, rs, rc,
     >                 d2a, dalpha, dbeta,
     >                 dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
          dra(il) = dr
c          
          nu = il*il
          ik = il+1
c          
          bx = ddlf(ik) * dr
          by = 0.0d0
          bz = -dlf(ik) * dr
     >       * dble(ik)
c
          dx = dx + bx * bc(nu) 
          dy = dy + by * bc(nu) 
          dz = dz + bz * bc(nu) 
        enddo
c
c
        do im = 1,shdeg
          call mk_lf_dlf(im, shdeg, rs, rc,
     >                   d2a, dalpha, dbeta,
     >                   dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            nu = (il*il + 2*(im-1)) + 1
            ik = il-im+1
            dr = dra(il)

            dw = ddlf(ik) * dr
            bx   = dw*dc
            bxp1 = dw*ds
            
            dw = (dlf(ik)/rs) * dr
     >         * dble(im) 
            by   =  dw*ds
            byp1 = -dw*dc
            
            dw = dlf(ik) * dr
     >         * dble(il+1)
            bz   = -dw*dc
            bzp1 = -dw*ds

            dx = dx + bx * bc(nu) + bxp1 * bc(nu+1) 
            dy = dy + by * bc(nu) + byp1 * bc(nu+1)   
            dz = dz + bz * bc(nu) + bzp1 * bc(nu+1) 
          enddo
        enddo
c
        deallocate(dalpha, dbeta)
        deallocate(dra, dlf, ddlf)
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

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	function XYZsph_bi0_fun
c		
c       Computes the dot product of 'be' and 'bc', avoiding the need
c       to store the entire contents of the 'be' coefficient array.
c             
c       input:
c          shdeg  INTEGER     max SH degree value
c          nb     INTEGER     number of coefficients
c          d2a    REAL*8      pre-computed array for mk_lf_dlf()
c          bc     REAL*8      coefficient array
c          p1     REAL*8      co-latitude
c          p2     REAL*8      longitude
c          ra     REAL*8      radius
c
c       output:
c          YZsph_bi0_fun  REAL*8
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        real*8 function XYZsph_bi0_fun(shdeg, nb,
     >                                 d2a, bc,
     >                                 p1, p2, ra,
     >                                 bex, bey, bez)
c
        implicit none
c
        integer shdeg, nb
        real*8 d2a(0:shdeg), bc(nb)
        real*8 p1, p2, ra
        real*8 bex, bey, bez 
c
        integer nu, il, im, ik
        real*8 rc, rs
        real*8 ds, dc, dr, dw
        real*8 bx, by, bz
        real*8 bxp1, byp1, bzp1
c
        real*8, allocatable :: dalpha(:), dbeta(:)
        real*8, allocatable :: dra(:), dlf(:), ddlf(:)
c 
c
#if defined(OMP_OFFLOAD_CPTP)
!$omp declare target
#endif
c
c
        allocate(dalpha(2:shdeg-1), dbeta(2:shdeg-1))
        allocate(dra(shdeg), dlf(shdeg+1), ddlf(shdeg+1)) 
c
        XYZsph_bi0_fun = 0.0d0 
c        
        rc = dcos(p1)
        rs = dsin(p1)
c
        im = 0
        call mk_lf_dlf(im, shdeg, rs, rc,
     >                 d2a, dalpha, dbeta,
     >                 dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
          dra(il) = dr
c          
          nu = il*il
          ik = il+1
c          
          bx = ddlf(ik) * dr
          by = 0.0d0
          bz = -dlf(ik) * dr
     >       * dble(ik)
c
          XYZsph_bi0_fun = XYZsph_bi0_fun
     >                   + (bex*bx+bey*by+bez*bz) * bc(nu) 
        enddo
c
c
        do im = 1,shdeg
          call mk_lf_dlf(im, shdeg, rs, rc,
     >                   d2a, dalpha, dbeta,
     >                   dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            nu = (il*il + 2*(im-1)) + 1
            ik = il-im+1
            dr = dra(il)
c
            dw = ddlf(ik) * dr
            bx   = dw*dc
            bxp1 = dw*ds
c            
            dw = (dlf(ik)/rs) * dr
     >         * dble(im) 
            by   =  dw*ds
            byp1 = -dw*dc
c            
            dw = dlf(ik) * dr
     >         * dble(il+1)
            bz   = -dw*dc
            bzp1 = -dw*ds
c
            XYZsph_bi0_fun = XYZsph_bi0_fun
     >                     + (bex*bx+bey*by+bez*bz) * bc(nu) 
     >                     + (bex*bxp1+bey*byp1+bez*bzp1) * bc(nu+1)  
          enddo
        enddo
c
        deallocate(dalpha, dbeta)
        deallocate(dra, dlf, ddlf)
c
        return
        end function XYZsph_bi0_fun


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine XYZsph_bi0_sub
c		
c       Computes the gradient of the weighted sum of squares (gj) and
c       the diagonal of the Hessian (hj). 
c
c       input:
c          shdeg    INTEGER     max SH degree value
c          nb       INTEGER     number of coefficients
c          d2a      REAL*8      pre-computed array for mk_lf_dlf()
c          p1       REAL*8      co-latitude
c          p2       REAL*8      longitude
c          ra       REAL*8      radius
c          dw_gj    REAL*8      multiplier for gj terms
c          dw_hj    REAL*8      multiplier for hj terms
c
c       output:
c          gj       REAL*8      gradient of the weighted sum of squares (nb)
c          hj       REAL*8      diagonal of the Hessian (nb)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine XYZsph_bi0_sub(shdeg, nb, d2a,
     >                            p1, p2, ra,
     >                            bex, bey, bez,
     >                            dw_gj, dw_hj,
     >                            gj, hj)
c
        implicit none
c
        integer shdeg, nb
        real*8 d2a(0:shdeg)
        real*8 p1, p2, ra
        real*8 bex, bey, bez 
        real*8 dw_gj, dw_hj
        real*8 gj(nb), hj(nb) 
c
        integer nu, il, im, ik
        real*8 rc, rs
        real*8 ds, dc, dr, dw
        real*8 bx, by, bz
        real*8 bxp1, byp1, bzp1
        real*8 be, bep1
c
        real*8, allocatable :: dalpha(:), dbeta(:)
        real*8, allocatable :: dra(:), dlf(:), ddlf(:)
c 
c
#if defined(OMP_OFFLOAD_SSQGH)
!$omp declare target
#endif
c
c
        allocate(dalpha(2:shdeg-1), dbeta(2:shdeg-1))
        allocate(dra(shdeg), dlf(shdeg+1), ddlf(shdeg+1))
c
        rc = dcos(p1)
        rs = dsin(p1)
c
        im = 0
        nu = 1
        call mk_lf_dlf(im, shdeg, rs, rc,
     >                 d2a, dalpha, dbeta,
     >                 dlf, ddlf)
        do il = 1,shdeg
          dr = ra**(il+2)
          dra(il) = dr
c          
          ik = il+1
c          
          bx = ddlf(ik) * dr
          by = 0.0d0
          bz = -dlf(ik) * dr
     >       * dble(ik)
c
          be = (bex*bx + bey*by + bez*bz)
c
!$OMP ATOMIC
          gj(nu) = gj(nu) + dw_gj*be
!$OMP ATOMIC
          hj(nu) = hj(nu) + dw_hj*be*be
c
          nu = nu+1
        enddo
c
c
        do im = 1,shdeg
          call mk_lf_dlf(im, shdeg, rs, rc,
     >                   d2a, dalpha, dbeta,
     >                   dlf, ddlf)
          dc = dcos(im*p2)
          ds = dsin(im*p2)
          do il = im,shdeg
            ik = il-im+1
            dr = dra(il)
c
            dw = ddlf(ik) * dr
            bx   = dw*dc
            bxp1 = dw*ds
c            
            dw = (dlf(ik)/rs) * dr
     >         * dble(im) 
            by   =  dw*ds
            byp1 = -dw*dc
c            
            dw = dlf(ik) * dr
     >         * dble(il+1)
            bz   = -dw*dc
            bzp1 = -dw*ds
c
            be = bex*bx + bey*by + bez*bz
            bep1 = bex*bxp1 + bey*byp1 + bez*bzp1 
c
!$OMP ATOMIC
            gj(nu) = gj(nu) + dw_gj*be
!$OMP ATOMIC
            gj(nu+1) = gj(nu+1) + dw_gj*bep1
!$OMP ATOMIC
            hj(nu) = hj(nu) + dw_hj*be*be
!$OMP ATOMIC
            hj(nu+1) = hj(nu+1) + dw_hj*bep1*bep1
c
            nu = nu+2
          enddo
        enddo
c
        deallocate(dalpha, dbeta)
        deallocate(dra, dlf, ddlf)
c
        return
        end subroutine XYZsph_bi0_sub