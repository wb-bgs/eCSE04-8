cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     module sph_wmam
c                       V. Lesur 16/09/2006
c
c     This module set the variable and define the model for the
c     anaomaly field data
c
c     All data are scalar
c
c  Number of parameters as organised on BE:
c  ilg*(ilg+2)                        internal lithosphere models
c
c    Global data:
c       ryg : real*8  : reference year for the model
c       ilg : integer : internal SH degree value (see below)
c
c    Subroutines:
c      init_sph_wmam.f    initialise the global variables
c      sub_sph_wmam.f     model subroutine
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        module sph_wmam
c
        implicit none
c
        integer ilg
        real*8 ryg,rag
c
        contains
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine init_sph_wmam
c
c       V. Lesur 16/09/2006
c
c  Set ilg, ryg (based on 1990) and rag
c
c  input:
c       shdeg : integer  : spherical harmonic degree
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine init_sph_wmam(shdeg)
c
        implicit none
c
        integer shdeg
c
        call dy2mjd(1990.d0,ryg)
        rag=6371.2d0
        ilg=shdeg
c
        return
        end subroutine init_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sub_sph_wmam
c
c       V. Lesur  16/09/2006
c
c   That is for a full non-lineair inversion of the lithosphere field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        recursive subroutine sub_sph_wmam(iof,nub,nb,bc,bp,be)
c
        implicit none
c
        character iof
        integer nub,nb
        integer i
        real*8 bp(*),be(*),bc(*)
        real*8, allocatable :: dwx(:),dwy(:),dwz(:)
c
        integer il,im
        real*8 dx,dy,dz,dd
        real*8 dxc,dyc,dzc
c
c  Damping
c
        if(nub.eq.100)then
c  setting
          allocate(dwx(1:nb),dwy(1:nb),dwz(1:nb))
          dwx(1:nb)=0.0d0
          dwy(1:nb)=0.0d0
          dwz(1:nb)=0.0d0
          be(1:nb)=0.0d0
c
c    calculate internal field component
          call XYZsph_bi0(ilg,rag,bp,dwx,dwy,dwz)
c
c    calculate field (lithos)
          dx=dot_product(dwx(1:nb),bc(1:nb))
          dy=dot_product(dwy(1:nb),bc(1:nb))
          dz=dot_product(dwz(1:nb),bc(1:nb))
c
c    calculate field unit vector (lithos)
          dd=dsqrt(bp(5)**2+bp(6)**2+bp(7)**2)
          dxc=bp(5)/dd
          dyc=bp(6)/dd
          dzc=bp(7)/dd
c
          dd=(dy*dzc-dyc*dz)**2
          dd=dd+(dx*dzc-dxc*dz)**2
          dd=dd+(dx*dyc-dxc*dy)**2
          dd=dsqrt(dd)
c
          if(iof.ne.'f')then
c    Calculate be(:)
            do i=1,nb
              be(i)=dwx(i)*((dx*dzc-dxc*dz)*dzc+(dx*dyc-dxc*dy)*dyc)
              be(i)=be(i)+dwy(i)
     >                    *((dy*dzc-dyc*dz)*dzc+(dy*dxc-dyc*dx)*dxc)
              be(i)=be(i)+dwz(i)
     >                    *((dz*dyc-dzc*dy)*dyc+(dz*dxc-dzc*dx)*dxc)
              be(i)=be(i)/dd
            enddo
          else
c    Full non-linear forward
            be(1:1)=dd
          endif
          deallocate(dwx,dwy,dwz)
c
c  Regular data
        else
c  setting
          allocate(dwx(1:nb),dwy(1:nb),dwz(1:nb))
          dwx(1:nb)=0.0d0
          dwy(1:nb)=0.0d0
          dwz(1:nb)=0.0d0
          be(1:nb)=0.0d0
c
c    calculate internal field component
          call XYZsph_bi0(ilg,rag,bp,dwx,dwy,dwz)
c
c    calulate full field (lithos + core)
          dx=bp(5)+dot_product(dwx(1:nb),bc(1:nb))
          dy=bp(6)+dot_product(dwy(1:nb),bc(1:nb))
          dz=bp(7)+dot_product(dwz(1:nb),bc(1:nb))
          dd=dsqrt(dx**2+dy**2+dz**2)
c
          if(iof.ne.'f')then
c    Calculate be(:)
            do i=1,nb
              be(i)=dwx(i)*dx
              be(i)=be(i)+dwy(i)*dy
              be(i)=be(i)+dwz(i)*dz
              be(i)=be(i)/dd
            enddo
          else
c    Full non-linear forward
            dx=dsqrt(bp(5)**2+bp(6)**2+bp(7)**2)
            be(1:1)=dd-dx
          endif
c
          deallocate(dwx,dwy,dwz)
        endif
c
        return
        end subroutine  sub_sph_wmam
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sub_sph_wmam_l
c
c       V. Lesur  16/09/2006
c
c   That is for a linearized inversion of the lithosphere field
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        recursive subroutine sub_sph_wmam_l(iof,nub,nb,bc,bp,be)
c
        implicit none
c
        character iof
        integer nub,nb
        integer i
        real*8 bp(*),be(*),bc(*)
        real*8, allocatable :: dwx(:),dwy(:),dwz(:)
c
        integer il,im
        real*8 dx,dy,dz,dd
        real*8 dxc,dyc,dzc
c
c  Damping
c
        if(nub.eq.100)then
c  setting
          allocate(dwx(1:nb),dwy(1:nb),dwz(1:nb))
          dwx(1:nb)=0.0d0
          dwy(1:nb)=0.0d0
          dwz(1:nb)=0.0d0
          be(1:nb)=0.0d0
c
c    calculate internal field component
          call XYZsph_bi0(ilg,rag,bp,dwx,dwy,dwz)
c       
c    calulate field (lithos & core)
          dx=dot_product(dwx(1:nb),bc(1:nb))
          dy=dot_product(dwy(1:nb),bc(1:nb))
          dz=dot_product(dwz(1:nb),bc(1:nb))
c
          dd=dsqrt(bp(5)**2+bp(6)**2+bp(7)**2)
          dxc=bp(5)/dd
          dyc=bp(6)/dd
          dzc=bp(7)/dd
c
          dd=(dy*dzc-dyc*dz)**2
          dd=dd+(dx*dzc-dxc*dz)**2
          dd=dd+(dx*dyc-dxc*dy)**2
          dd=dsqrt(dd)
c
c         if(iof.ne.'f')then
c    Calculate be(:)
            do i=1,nb
              be(i)=dwx(i)*((dx*dzc-dxc*dz)*dzc+(dx*dyc-dxc*dy)*dyc)
              be(i)=be(i)+dwy(i)
     >                    *((dy*dzc-dyc*dz)*dzc+(dy*dxc-dyc*dx)*dxc)
              be(i)=be(i)+dwz(i)
     >                    *((dz*dyc-dzc*dy)*dyc+(dz*dxc-dzc*dx)*dxc)
              be(i)=be(i)/dd
            enddo
c         else
c    Full non-linear forward
c           be(1:1)=dd
c         endif
          if(iof.eq.'f')be(1:nb)=be(1:nb)*bc(1:nb)
          deallocate(dwx,dwy,dwz)
c
c  Regular data
        else
c
c    calculate full field (core)
          dd=dsqrt(bp(5)**2+bp(6)**2+bp(7)**2)
c
c  setting
          allocate(dwx(1:nb),dwy(1:nb),dwz(1:nb))
          dwx(1:nb)=0.0d0
          dwy(1:nb)=0.0d0
          dwz(1:nb)=0.0d0
          be(1:nb)=0.0d0
c
c    calculate internal field component
          call XYZsph_bi0(ilg,rag,bp,dwx,dwy,dwz)
c
c    Calculate be(:)
          do i=1,nb
            be(i)=dwx(i)*bp(5)
            be(i)=be(i)+dwy(i)*bp(6)
            be(i)=be(i)+dwz(i)*bp(7)
            be(i)=be(i)/dd
          enddo
          if(iof.eq.'f')be(1:nb)=be(1:nb)*bc(1:nb)
c
          deallocate(dwx,dwy,dwz)
        endif
c
        return
        end subroutine  sub_sph_wmam_l
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function wmam_norm(i,nub,mv)
c
         implicit none
c
         integer i,nub(*)
         real*8 mv(*)
c
         real*8 dc,ae
c
         real*8 lesur_norm,l2_norm,l1_norm
         external lesur_norm,l2_norm,l1_norm
c
         select case (nub(i))
         case (1)
           dc=0.6d0
           ae=0.5d0
           wmam_norm=l2_norm(i,nub,mv)
c          wmam_norm=lesur_norm(dc,ae,i,mv)
         case default
           wmam_norm=l2_norm(i,nub,mv)
c          wmam_norm=l1_norm(i,nub,mv)
         end select
c
         return
         end function wmam_norm
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function wmam_var(nub,npt,np,var,mv,wgh)
c
         implicit none
c
         integer npt,np,nub(*)
         real*8 mv(*),var,wgh(*)
c
         integer i,j,k
         real*8 vv,dv
         real*8 dc,ae
c
         real*8 lesur_var2,l2_var,l1_var
         external lesur_var2,l2_var,l1_var
c
         if(np.ge.0)then
           vv=var
           do i=1,np
             dv=vv
             j=0
             k=1
c
             select case (nub(i))
             case (1)
               dc=0.6d0
               ae=0.5d0
               vv=l2_var(nub(i),j,k,dv,mv(i),wgh(i))
c              vv=lesur_var2(dc,ae,j,k,dv,mv(i),wgh(i))
             case default
               vv=l2_var(nub(i),j,k,dv,mv(i),wgh(i))
c              vv=l1_var(nub(i),j,k,dv,mv(i),wgh(i))
             end select
c
           enddo
         else
          vv=0.0d0
          do i=1,-np
            vv=vv+mv(i)
          enddo
         endif
c
         wmam_var=vv
c
         return
         end function wmam_var
ccccccccccccccccccccccccccccccccccc
        end module sph_wmam
