cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Subroutine XYZsph_bi0
c		Vincent Lesur 21/01/2009
c
c    modified:
c       07.01.2010  V.Lesur :
c         1- reference radius ra as parameter
c
c       For a all spherical harmonique between degree 1 and nlie
c       and a  colatitude, longitude, radius computes for
c       internal sources: 
c       X component basis  (X horizontal North)
c       Y component basis  (Y horizontal East)
c       Z component basis  (Z Vertical Down) 
c
c       These basis is the derivative of a potential 
c       defined on the sphere of radius "ra"  
c       (See Foundations of Geomagnetism, Backus 1996 page 110 and 125)
c      
c       No time dependence
c
c       Spherical harmonics are set by increasing degree "l", for 
c       a given degree by increasing order "m" (l=0,m=0 excluded)
c       nu=(l**2+2*(m-1))   FOR cosine terms if m.ne.0
c       nu=(l**2+2*m-1)     For sine terms or m=0
c
c       colat,long,radius given in degree,degree,km
c
c       input:
c         ilg  INTEGER   max SH degree
c         nd   INTEGER   space dimension 
c         rag  REAL*8    reference radius
c         d2a  REAL*8    pre-computed d2a array for mklf_F2()
c         pos  REAL*8    positions in space and time (colat,long,radius)
c
c       output:
c         bx   REAL*8   X Base funtion value
c         by   REAL*8   Y Base funtion value
c         bz   REAL*8   Z Base funtion value
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine XYZsph_bi0(ilg, nd, rag,
     >                        d2a, pos,
     >                        bx, by, bz)
c
c       original unoptimised version 
c
        implicit none
c
        integer ilg, nd
        real*8 rag
        real*8 d2a(0:ilg), pos(nd+1)
        real*8 bx(*),by(*),bz(*)
c
        integer nu,il,im,ik
        real*8 rc,rs,dw,d2r
        real*8 ds,dc,ra_div_pos3
        real*8 dr,p1,p2,p3
c
        real*8, allocatable :: dlf(:),ddlf(:)
        real*8, allocatable :: dra(:)
c
c
        allocate(dlf(ilg+1),ddlf(ilg+1),dra(ilg))
c
        d2r=4.d0*datan(1.d0)/180.d0
c
        p1=pos(1)
        p2=pos(2)
        p3=pos(3)

        rc = dcos(p1*d2r)
        rs = dsin(p1*d2r)

        ra_div_pos3 = rag / p3
c
c   im=0
        im=0
        call mk_lf_dlf(im,ilg,rs,rc,d2a,dlf,ddlf)
        do il=1,ilg
          dr = ra_div_pos3**(il+2)
          dra(il) = dr
c          
          nu = il*il
          ik = il+1
c          
          bx(nu) = ddlf(ik) * dr
          by(nu) = 0.0d0
          bz(nu) = -dlf(ik) * dr
     >           * dble(ik)
        enddo

c   im.ne.0
        do im=1,ilg
          call mk_lf_dlf(im,ilg,rs,rc,d2a,dlf,ddlf)
          dc = dcos(im*p2*d2r)
          ds = dsin(im*p2*d2r)
          do il=im,ilg
            nu = (il*il + 2*(im-1)) + 1
            ik = il-im+1
            dr = dra(il)

            dw = ddlf(ik) * dr
            bx(nu)   = dw*dc
            bx(nu+1) = dw*ds
            
            dw = (dlf(ik)/rs) * dr
     >         * dble(im) 
            by(nu)   =  dw*ds
            by(nu+1) = -dw*dc
            
            dw = dlf(ik) * dr
     >         * dble(il+1)
            bz(nu)   = -dw*dc
            bz(nu+1) = -dw*ds
          enddo
        enddo
c
        deallocate(dlf,ddlf,dra)
c
        return
        end

        subroutine XYZsph_bi0_fun(ilg, nb, nd, rag,
     >                            d2a, bc, pos, bedotbc,
     >                            bex, bey, bez)
c
c       this variant takes 'bex', 'bey', and 'bez' as inputs and 
c       computes 'bedotbc' as the dot product of 'be' and 'bc', but
c       does not store the values of 'be' 
c
        implicit none
c
        integer ilg, nb, nd
        real*8 rag
        real*8 d2a(0:ilg), bc(nb), pos(nd+1)
        real*8 bedotbc
        real*8 bex, bey, bez 
c
        integer nu,il,im,ik
        real*8 rc,rs,dw,d2r
        real*8 ds,dc,ra_div_pos3
        real*8 dr,p1,p2,p3
        real*8 bx,by,bz,bxp1,byp1,bzp1
c
        real*8, allocatable :: dlf(:),ddlf(:)
        real*8, allocatable :: dra(:)
c 
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
        allocate(dlf(ilg+1),ddlf(ilg+1),dra(ilg))

        bedotbc = 0.0d0 
        d2r=4.d0*datan(1.d0)/180.d0

        p1=pos(1)
        p2=pos(2)
        p3=pos(3)

        rc = dcos(p1*d2r)
        rs = dsin(p1*d2r)

        ra_div_pos3 = rag / p3
c
c   im=0
        im=0
        call mk_lf_dlf(im,ilg,rs,rc,d2a,dlf,ddlf)
        do il=1,ilg
          dr = ra_div_pos3**(il+2)
          dra(il) = dr
c          
          nu = il*il
          ik = il+1
c          
          bx = ddlf(ik) * dr
          by = 0.0d0
          bz = -dlf(ik) * dr
     >           * dble(ik)
          bedotbc = bedotbc + (bex*bx+bey*by+bez*bz) *bc(nu) 

        enddo

c   im.ne.0
        do im=1,ilg
          call mk_lf_dlf(im,ilg,rs,rc,d2a,dlf,ddlf)
          dc = dcos(im*p2*d2r)
          ds = dsin(im*p2*d2r)
          do il=im,ilg
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

            bedotbc = bedotbc + (bex*bx+bey*by+bez*bz) * bc(nu) 
     >                + (bex*bxp1+bey*byp1+bez*bzp1) * bc(nu+1)  
          enddo
        enddo
c
        deallocate(dlf,ddlf,dra)
c
        return
        end

        subroutine XYZsph_bi0_sub(ilg, nb, nd, rag,
     >                            d2a, be, pos,
     >                            bex, bey, bez)
c
c       this variant takes 'bex', 'bey', and 'bez' as inputs and 
c       computes 'be' 
c
        implicit none
c
        integer ilg, nb, nd
        real*8 rag
        real*8 d2a(0:ilg), be(nb), pos(nd+1)
        real*8 bex, bey, bez 
c
        integer nu,il,im,ik
        real*8 rc,rs,dw,d2r
        real*8 ds,dc,ra_div_pos3
        real*8 dr,p1,p2,p3
        real*8 bx,by,bz,bxp1,byp1,bzp1
c
        real*8, allocatable :: dlf(:),ddlf(:)
        real*8, allocatable :: dra(:)
c 
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
        allocate(dlf(ilg+1),ddlf(ilg+1),dra(ilg))

        d2r=4.d0*datan(1.d0)/180.d0
c
        p1=pos(1)
        p2=pos(2)
        p3=pos(3)

        rc = dcos(p1*d2r)
        rs = dsin(p1*d2r)

        ra_div_pos3 = rag / p3
c
c   im=0
        im=0
        call mk_lf_dlf(im,ilg,rs,rc,d2a,dlf,ddlf)
        do il=1,ilg
          dr = ra_div_pos3**(il+2)
          dra(il) = dr
c          
          nu = il*il
          ik = il+1
c          
          bx = ddlf(ik) * dr
          by = 0.0d0
          bz = -dlf(ik) * dr
     >           * dble(ik)
          be(nu) = bex*bx+bey*by+bez*bz 

        enddo

c   im.ne.0
        do im=1,ilg
          call mk_lf_dlf(im,ilg,rs,rc,d2a,dlf,ddlf)
          dc = dcos(im*p2*d2r)
          ds = dsin(im*p2*d2r)
          do il=im,ilg
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

            be(nu) = bex*bx+bey*by+bez*bz 
            be(nu+1) = bex*bxp1+bey*byp1+bez*bzp1 

          enddo
        enddo
c
        deallocate(dlf,ddlf,dra)   
c
        return
        end

        subroutine XYZsph_bi0_sample(ilg, nb, nd, rag,
     >                               d2a, bc, pos,
     >                               dx, dy, dz)
c
c       this variant is only called for sample points 
c       before calling XYZsph_bi0_fun or XYZsph_bi0_sub 
c       it computes dx = bx.bc, dy = by.bc and dz = bz.bc 
c
        implicit none
c
        integer ilg, nb, nd
        real*8 rag
        real*8 d2a(0:ilg), bc(nb), pos(nd+1)
        real*8 dx, dy, dz 
c
        integer nu,il,im,ik
        real*8 rc,rs,dw,d2r
        real*8 ds,dc,ra_div_pos3
        real*8 dr,p1,p2,p3
        real*8 bx,by,bz,bxp1,byp1,bzp1
c
        real*8, allocatable :: dlf(:),ddlf(:)
        real*8, allocatable :: dra(:) 
c 
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
        allocate(dlf(ilg+1),ddlf(ilg+1),dra(ilg))
c
        d2r=4.d0*datan(1.d0)/180.d0

        dx = 0.0d0
        dy = 0.0d0 
        dz = 0.0d0 

        p1=pos(1)
        p2=pos(2)
        p3=pos(3)

        rc = dcos(p1*d2r)
        rs = dsin(p1*d2r)

        ra_div_pos3 = rag / p3
c
c   im=0
        im=0
        call mk_lf_dlf(im,ilg,rs,rc,d2a,dlf,ddlf)
        do il=1,ilg
          dr = ra_div_pos3**(il+2)
          dra(il) = dr
c          
          nu = il*il
          ik = il+1
c          
          bx = ddlf(ik) * dr
          by = 0.0d0
          bz = -dlf(ik) * dr
     >           * dble(ik)

          dx = dx + bx * bc(nu) 
          dy = dy + by * bc(nu) 
          dz = dz + bz * bc(nu) 
        enddo

c   im.ne.0
        do im=1,ilg
          call mk_lf_dlf(im,ilg,rs,rc,d2a,dlf,ddlf)
          dc = dcos(im*p2*d2r)
          ds = dsin(im*p2*d2r)
          do il=im,ilg
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
        deallocate(dlf,ddlf,dra)
c
        return
        end
