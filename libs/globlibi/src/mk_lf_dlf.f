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
c         nm/nl         order and degree max
c         dc/ds         cos(colatitude)/sin(colatitude)
c       output:
c         dlf           legendre function from nm to nl
c         ddlf          derivative of legendre function from nm to nl
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine mk_lf_dlf(nm,nl,ds,dc,dlf,ddlf)
c
        implicit none
c
        integer nm,nl,il,jl
        real*8 ds,dc,dlf(*),ddlf(*),d1,d2
c
        call mklf_F(nm,nl,ds,dc,dlf)
c
        if (ds.eq.0.0d0) then
          if(nm.ne.1) then 
            do il=1,nl+1
              ddlf(il)=0.0d0
            enddo
          else
            ddlf(1)=-1.0d0
            ddlf(2)=-dsqrt(3.0d0) 
            do il=3,nl+1
              d1=(2*il-1)/dsqrt(dble(il*il-1))
              d2=dsqrt(dble((il-1)**2-1)/(il*il-1))
              ddlf(il)=d1*ddlf(il-1)-d2*ddlf(il-2)
            enddo
          endif
        else
          do il=nl,nm,-1
            jl=il-nm+1
            if(il-nm.gt.0)then
              d1=dsqrt(dble((il-nm)*(il+nm)))
              d2=dble(il) 
              ddlf(jl)=(d2*dc*dlf(jl)-d1*dlf(jl-1))/ds
            else
              ddlf(1)=nm*dc*dlf(1)/ds
            endif
          enddo
        endif
c
        return
        end
