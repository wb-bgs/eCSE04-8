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
c         nm/ilg        order and degree max
c         dc/ds         cos(colatitude)/sin(colatitude)
c         d2a           pre-computed array for mklf_F2()
c       output:
c         dlf           legendre function from nm to nl
c         ddlf          derivative of legendre function from nm to nl
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mk_lf_dlf(nm, ilg, ds, dc, d2a, dlf, ddlf)
c
        implicit none
c
        integer nm, ilg
        real*8 ds, dc, d2a(0:ilg)
        real*8 dlf(ilg+1), ddlf(ilg+1)
c
        integer il, jl
        real*8 d1, d2
c
c
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
        call mklf_F2(nm, ilg, ds, dc, d2a, dlf)
c
        if (ds.eq.0.0d0) then

            if (nm.ne.1) then 
                ddlf(1:ilg+1) = 0.0d0
            else
                ddlf(1) = -1.0d0
                ddlf(2) = -dsqrt(3.0d0) 
                do il=3,ilg+1
                    jl = il*il-1
                    d1 = (2*il-1)/dsqrt(dble(jl))
                    d2 = dsqrt(dble((il-1)**2-1)/jl)
                    ddlf(il) = d1*ddlf(il-1) - d2*ddlf(il-2)
                enddo
            endif

        else

            do il=ilg,nm+1,-1
                jl = il - nm + 1
                d1 = dsqrt(dble((il-nm)*(il+nm)))
                d2 = dble(il) 
                ddlf(jl) = (d2*dc*dlf(jl) - d1*dlf(jl-1)) / ds
            enddo
            ddlf(1) = nm * dc * dlf(1)/ds

        endif
c
        return
        end
