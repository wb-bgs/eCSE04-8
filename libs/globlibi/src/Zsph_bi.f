cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Subroutine Zsph_bi
c		Vincent Lesur 30/04/2002
c
c       For a all spherical harmonique between nlis,nlie
c       and a  colatitude, longitude, radius and time compute
c       the value of the Z component basis  (Z vertical DOWN) for
c       internal sources only. This basis 
c       is the vertical derivative of an potential defined on the sphere 
c       of radius "ra"
c       (See Foundations of Geomagnetism, Backus 1996 page 110 and 125)
c       Time dependence of the potential is introduce via a polynomial in
c       time.
c
c       Spherical harmonics are set by increasing degree "l", for 
c       a given degree by increasing order "m" (l=0,m=0 excluded) and for
c       each lm paires by increasing polynomial degree in time: 
c       nu=(l**2+2*(m-1))*nti+it   FOR cosine terms if m.ne.0
c       nu=(l**2+2*m-1)*nti+it     For sine terms or m=0
c
c       colat,long,radius given in degree,degree,km
c
c       input:
c         ry            reference year
c         nlis          start of degree internal
c         nlie          end of degree internal
c         nti           number of time dependent parameters internal
c         pos           positions in space and time (colat,long,radius,jd,refjd)
c
c       output:
c         be		Base funtion value
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine Zsph_bi(ry,nlis,nlie,nti,pos,be)
c
        implicit none
c
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
        integer nu,nti,il,im,nlis,nlie,nus,it
        real*8 ry,pos(*),be(*),rc,rs,ra,dw
        real*8 dc,ds,dt,p1,p2,p4
        real*8, allocatable :: dlf(:)
c
        allocate(dlf(1:nlie+1))
c
        ra=6371.2d0
c
        p1=pos(1)*D2R
        p2=pos(2)*D2R
        p4=pos(4)

        dt=(p4-ry)/365.25d0
        rc=dcos(p1)
        rs=dsin(p1)
        nus=(nlis*nlis-1)*nti
c
c Internal
c   im=0
        im=0
        call mklf_F(im,nlie,rs,rc,dlf)
        do il=nlis,nlie
          nu=(il*il-1)*nti-nus
          be(nu+1)=-dlf(il-im+1)*(il+1)*(ra/pos(3))**(il+2)
          do it=2,nti
            be(nu+it)=be(nu+1)*dt**(it-1)
          enddo
        enddo
c   im.ne.0
        do im=1,nlie
          call mklf_F(im,nlie,rs,rc,dlf)
          dc=dcos(im*p2)
          ds=dsin(im*p2)
          do il=max0(im,nlis),nlie
            dw=dlf(il-im+1)*(il+1)*(ra/pos(3))**(il+2)
            nu=(il*il+2*(im-1))*nti-nus
            be(nu+1)=-dw*dc
            be(nu+nti+1)=-dw*ds
            do it=2,nti
              dw=dt**(it-1)
              be(nu+it)=be(nu+1)*dw
              be(nu+nti+it)=be(nu+nti+1)*dw
            enddo
          enddo
        enddo
c
        deallocate(dlf)
        return
        end
