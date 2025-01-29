cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Subroutine Xsph_bi
c		Vincent Lesur 30/04/2002
c
c       For a all spherical harmonique between degree nlis and nlie
c       and a  colatitude, longitude, radius and time computes 
c       the value of the X component basis  (X horizontal North) for 
c       internal sources. 
c       This basis is the derivative along the (-colatitude) of a potential 
c       defined on the sphere of radius "ra"  
c       (See Foundations of Geomagnetism, Backus 1996 page 110 and 125)
c       Time dependence of the potential is introduced via a polynomial in
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
c         nd            space dimension
c         nb            number of base element 
c         nlie          end of degree internal
c         pos           positions in space and time (colat,long,radius,jd,refjd)
c
c       output:
c         be		Base funtion value
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine Xsph_bi(nd, nb, nlie, pos, be)
c
        use kernels, only : RAG, D2R
c
        implicit none
c
        integer nd, nb, nlie
        real*8 pos(1:nd+1), be(1:nb)
c
        integer nlis, nti
        integer nu, il, im, nus, it
        real*8 rc, rs, dw
        real*8 ds, dc, dt
        real*8 p1, p2, p3
        real*8, allocatable :: ddlf(:)
c
c
        nlis = 1
        nti = 1
        allocate(ddlf(1:nlie+1))
c
        p1 = pos(1)*D2R
        p2 = pos(2)*D2R
        p3 = RAG/pos(3)
c
        dt = pos(4)/365.25d0
        rc = dcos(p1)
        rs = dsin(p1)
        nus = (nlis*nlis - 1)*nti
c
        im=0
        call mkdlf_F(im, nlie, rs, rc, ddlf)
        do il = nlis,nlie
            nu = (il*il-1)*nti - nus
            be(nu+1) = ddlf(il-im+1)*(p3**(il+2))
            do it = 2,nti
                be(nu+it) = be(nu+1)*(dt**(it-1))
            enddo
        enddo
c 
        do im = 1,nlie
            call mkdlf_F(im, nlie, rs, rc, ddlf)
            dc = dcos(im*p2)
            ds = dsin(im*p2)
            do il = max0(im,nlis),nlie
                dw = ddlf(il-im+1)*(p3**(il+2))
                nu = (il*il + 2*(im-1))*nti - nus
                be(nu+1) = dw*dc
                be(nu+nti+1) = dw*ds
                do it = 2,nti
                    dw = dt**(it-1)
                    be(nu+it) = be(nu+1)*dw
                    be(nu+nti+it) = be(nu+nti+1)*dw
                enddo
            enddo
        enddo
c
        deallocate(ddlf)
c
        return
        end