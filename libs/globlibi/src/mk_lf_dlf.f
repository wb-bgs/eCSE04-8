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
c
c       input:
c         ilg  INTEGER  degree max
c         nm   INTEGER  order
c         ds   REAL*8  sin(colatitude)
c         dc   REAL*8  cos(colatitude)
c         d2a, d3a, dalphaa, dbetaa  REAL*8
c           arrays holding pre-computed values for mklf_F2()
c         d4a, d5a, d6a  REAL*8
c           arrays holding pre-computed values for mk_lf_dlf()
c       output:
c         dlf   REAL*8  legendre function from nm to nl
c         ddlf  REAL*8  derivative of legendre function from nm to nl
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mk_lf_dlf(nm, ilg, d2a, d3a,
     >                       dalphaa, dbetaa,
     >                       d4a, d5a, d6a,
     >                       ds, dc, dlf, ddlf)
c
        implicit none
c
        integer nm, ilg
        real*8 d2a(0:ilg), d3a(0:ilg)
        real*8 dalphaa(*), dbetaa(*)
        real*8 d4a(1:ilg+1), d5a(3:ilg+1), d6a(*)
        real*8 ds, dc
        real*8 dlf(*), ddlf(*)
c
        integer il, jl
        integer, save :: k=0
        
c
        call mklf_F2(nm, ilg, d2a, d3a,
     >               dalphaa, dbetaa,
     >               ds, dc, dlf)

        if (ds .eq. 0.0d0) then

          if (nm .ne. 1) then 
            ddlf(1:ilg+1) = 0.0d0
          else
            ddlf(1) = d4a(1)
            ddlf(2) = d4a(2)
            do il=3,ilg+1
              ddlf(il) = d4a(il)*ddlf(il-1) - d5a(il)*ddlf(il-2)
            enddo
          endif

        else

          if (nm .eq. 0) then
            k = 1
          endif

          do il = ilg,nm,-1
            jl = il - nm + 1
            if (il-nm .gt. .0) then
              ddlf(jl) = (dble(il)*dc*dlf(jl) - d6a(k)*dlf(jl-1)) / ds
            else
              ddlf(1) = nm * dc * dlf(1)/ds
            endif
            k = k+1
          enddo

        endif
c
        return
        end