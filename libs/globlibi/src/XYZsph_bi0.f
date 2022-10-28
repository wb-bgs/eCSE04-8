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
c         rag  REAL*8    reference radius
c         pos  REAL*8    positions in space and time (colat,long,radius)
c         d2a  REAL*8    pre-computed d2a array for mklf_F2()
c
c       output:
c         bx   REAL*8   X Base funtion value
c         by   REAL*8   Y Base funtion value
c         bz   REAL*8   Z Base funtion value
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine XYZsph_bi0(ilg, rag, pos, d2a,
     >                        bx, by, bz)
c
        implicit none
c
        integer ilg
        real*8 rag
        real*8 pos(*), d2a(0:ilg)
        real*8 bx(*),by(*),bz(*)
        real*8, allocatable :: dlf(:),ddlf(:)
        real*8, allocatable :: dra(:)
c
        integer nu,il,im,ik
        real*8 rc,rs,dw
        real*8 ds,dc,dcosd,dsind,ra_div_pos3
        real*8 dr,p1,p2,p3
c
        external dcosd,dsind
        allocate(dlf(ilg+1),ddlf(ilg+1),dra(ilg))
c
        p1=pos(1)
        p2=pos(2)
        p3=pos(3)

        rc = dcosd(p1)
        rs = dsind(p1)

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
          dc = dcosd(im*p2)
          ds = dsind(im*p2)
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

        deallocate(dlf,ddlf,dra)
c
        return
        end