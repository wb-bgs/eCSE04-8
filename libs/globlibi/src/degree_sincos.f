ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   degree_sincos functions
c
c   Calulate de cosinus and sinus of an angle given in degree
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function dsind(theta)
c
         implicit none
c
         real*8 theta,d2r
c
         d2r=4.d0*datan(1.d0)/180.d0
         dsind=dsin(theta*d2r)
c
         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function dcosd(theta)
c
         implicit none
c
         real*8 theta,d2r
c
         d2r=4.d0*datan(1.d0)/180.d0
         dcosd=dcos(theta*d2r)
c
         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function dtand(theta)
c
         implicit none
c
         real*8 theta,d2r
c
         d2r=4.d0*datan(1.d0)/180.d0
         dtand=dtan(theta*d2r)
c
         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function dasind(val)
c
         implicit none
c
         real*8 val,d2r
c
         d2r=4.d0*datan(1.d0)/180.d0
         dasind=dasin(val)/d2r
c
         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function dacosd(val)
c
         implicit none
c
         real*8 val,d2r
c
         d2r=4.d0*datan(1.d0)/180.d0
         dacosd=dacos(val)/d2r
c
         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function datand(val)
c
         implicit none
c
         real*8 val,d2r
c
         d2r=4.d0*datan(1.d0)/180.d0
         datand=datan(val)/d2r
c
         return
         end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         real*8 function datan2d(val1,val2)
c
         implicit none
c
         real*8 val1,val2,d2r
c
         d2r=4.d0*datan(1.d0)/180.d0
         datan2d=datan2(val1,val2)/d2r
c
         return
         end
