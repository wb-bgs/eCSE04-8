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
c         nm/shdeg      order and degree max
c         dc/ds         cos(colatitude)/sin(colatitude)
c         d2a           pre-computed array for mklf_F2()
c       output:
c         dlf           legendre function from nm to nl
c         ddlf          derivative of legendre function from nm to nl
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mk_lf_dlf(nm, shdeg, ds, dc, d2a, dlf, ddlf)
c
        implicit none
c
        integer nm, shdeg
        real*8 ds, dc, d2a(0:shdeg)
        real*8 dlf(shdeg+1), ddlf(shdeg+1)
c
        integer il, jl
        real*8 d1, d2
c
c
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
c
        call mklf_F2(nm, shdeg, ds, dc, d2a, dlf)
c
        if (ds.eq.0.0d0) then
c
            if (nm.ne.1) then 
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
            do il=shdeg,nm+1,-1
                jl = il - nm + 1
                d1 = dsqrt(dble((il-nm)*(il+nm)))
                d2 = dble(il) 
                ddlf(jl) = (d2*dc*dlf(jl) - d1*dlf(jl-1)) / ds
            enddo
            ddlf(1) = nm * dc * dlf(1)/ds
c
        endif
c
        return
        end
