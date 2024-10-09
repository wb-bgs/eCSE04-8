cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c         d2a           pre-computed array for populating dlf
c       output:
c         dlf           legendre function from im to nl
c         ddlf          derivative of legendre function from im to nl
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mk_lf_dlf(im, shdeg, rs, rc, 
     >                       d2a, dalpha, dbeta,
     >                       dlf, ddlf)
c
        implicit none
c
        integer im, shdeg
        real*8 rs, rc
        real*8 d2a(0:shdeg)
        real*8 dalpha(2:shdeg-1), dbeta(2:shdeg-1)
        real*8 dlf(shdeg+1), ddlf(shdeg+1)
c
        integer d0, il, jl
        real*8 d1, d2
c
c
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
c
c  Initialise dlf array
        d1 = rs
        if (d1.ne.0.0d0) then
          d1 = d1**im
        else
          if (im.eq.0) d1 = 1.d0
        endif
c
        dlf(1) = d1*d2a(im)                         ! leg. func. (m,m)
        dlf(2) = dlf(1) * rc * dsqrt(dble(2*im+1))  ! leg. func. (m+1,m)
c
        do il=2,shdeg-im                     ! l=im+il-1
          d0 = il+2*im
          d1 = dble((il-1) * (d0-1))         ! (l-m)*(l+m)
          d2 = dble(il * d0)                 ! (l-m+1)*(l+m+1)	
          dbeta(il) = dsqrt(d1/d2)           ! recurrence coeff.
          d1 = dble(2*(il+im)-1)             ! 2l+1
          dalpha(il) = d1/dsqrt(d2)          !
        enddo
c
        do il=2,shdeg-im                     ! l=im+il-1
          ! leg. func. (im+il-1,im) 
          dlf(il+1) = dalpha(il)*dlf(il)*rc - dbeta(il)*dlf(il-1)
        enddo
c
c
c  Initialise ddlf array
        if (rs.eq.0.0d0) then
c
          if (im.ne.1) then 
            ddlf(1:shdeg+1) = 0.0d0
          else
            ddlf(1) = -1.0d0
            ddlf(2) = -dsqrt(3.0d0) 
            do il=3,shdeg+1
              jl = il*il-1
              d1 = (2*il-1)/dsqrt(dble(jl))
              d2 = dsqrt(dble((il-1)**2-1)/jl)
              ddlf(il) = d1*ddlf(il-1) - d2*ddlf(il-2)
            enddo
          endif
c
        else
c
          do il=shdeg,im+1,-1
            jl = il - im + 1
            d1 = dsqrt(dble((il-im)*(il+im)))
            d2 = dble(il) 
            ddlf(jl) = (d2*rc*dlf(jl) - d1*dlf(jl-1)) / rs
          enddo
          ddlf(1) = im * rc * dlf(1)/rs
c
        endif
c
        return
        end