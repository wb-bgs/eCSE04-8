ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sub_sph_wmam_l
c
c       V. Lesur  16/09/2006
c
c   That is for a linearized inversion of the lithosphere field.
c   This subroutine is called from either wmam_sub() or wmam_fun().
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine sub_sph_wmam_l(issampt, shdeg, nb, nb2, nd,
     >                            d2a, bc, bp, be,
     >                            bedotbc, fun)
c
	implicit none
c
        real*8, parameter :: RAG = 6371.2d0
c
	logical issampt  
	integer shdeg, nb, nb2, nd
	real*8 d2a(0:shdeg), bc(nb), bp(nd+1), be(nb2)
	real*8 bedotbc 
	logical fun
c
	real*8 dx, dy, dz, dd
c
 	real*8 dxbey, dxbez
	real*8 dybex, dybez
	real*8 dzbex, dzbey
c
	real*8 xy_c, xz_c
	real*8 yx_c, yz_c
	real*8 zx_c, zy_c
c
	real*8 bex, bey, bez
	real*8 bex2, bey2, bez2
c
c
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
c 
c  if 'fun' is set to .true. then this routine is being called from wmam_fun() 
c  and returns the dot product of 'be' and 'bc' in 'bedotbc', but 'be' does not
c  contain valid values
c
c  if 'fun' is set to .false. then this routine is being called from wmam_sub()
c  and returns the array 'be', and 'bedotbc' does not contain a valid value.   
c     
c  calculate internal field component
	if (issampt) then
c  this is a sampling point
          call XYZsph_bi0_sample(shdeg, nb, nd, RAG,
     >                           d2a, bc, bp,
     >                           dx, dy, dz)
	
	  bex = bp(5)
	  bey = bp(6)
	  bez = bp(7)
c
	  dxbey = dx*bey
	  dxbez = dx*bez
	  dybex = dy*bex
	  dybez = dy*bez
	  dzbex = dz*bex
	  dzbey = dz*bey
c
	  xy_c = dxbey - dybex
	  xz_c = dxbez - dzbex
	  yx_c = -xy_c
	  yz_c = dybez - dzbey
	  zx_c = -xz_c
	  zy_c = -yz_c
c
	  dd = dsqrt(yz_c**2 + xz_c**2 + xy_c**2)
c
	  bex2 = (xz_c*bez + xy_c*bey) / dd
	  bey2 = (yz_c*bez + yx_c*bex) / dd
	  bez2 = (zy_c*bey + zx_c*bex) / dd
	
	  bex = bex2
	  bey = bey2
	  bez = bez2
	
	else
c  this is a data point 
	
	  bex = bp(5)
	  bey = bp(6)
	  bez = bp(7)
	
	endif 
	
	if (fun) then 
c  computes 'bedotbc' as dot product of 'be' and 'bc' 
          call XYZsph_bi0_fun(shdeg, nb, nd, RAG,
     >                        d2a, bc, bp, bedotbc,
     >                        bex, bey, bez)
        else
c  computes 'be' 
          call XYZsph_bi0_sub(shdeg, nb, nd, RAG,
     >                        d2a, be, bp,
     >                        bex, bey, bez)
	end if 
c
	return
	end subroutine sub_sph_wmam_l
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine wmam_sub
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine wmam_sub(issampt, shdeg, nb, nd,
     >                      d2a, bc, bp, be)
c
	implicit none
c
        logical issampt
	integer shdeg, nb, nd
	real*8 d2a(0:shdeg), bc(nb), bp(nd+1), be(nb)
	real*8 dummybedotbc
c
        integer nb2
c
c
#ifdef OMP_OFFLOAD
!$omp declare target
#endif
c
c
        nb2 = nb
	call sub_sph_wmam_l(issampt, shdeg, nb, nb2, nd,
     >                      d2a, bc, bp, be,
     >                      dummybedotbc, .false.)
c 
	return
	end subroutine wmam_sub
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    function wmam_fun
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	real*8 function wmam_fun(issampt, shdeg, nb, nd,
     >                           d2a, bc, bp)
c
	implicit none
c
        logical issampt
	integer shdeg, nb, nd
	real*8 d2a(0:shdeg), bc(nb), bp(nd+1), bedotbc
	real*8 dummybe(1)
c
        integer nb2
c
c
#ifdef OMP_OFFLOAD 
!$omp declare target
#endif
c
c
        nb2=1
	call sub_sph_wmam_l(issampt, shdeg, nb, nb2, nd,
     >                      d2a, bc, bp, dummybe,
     >                      bedotbc, .true.)
c
	wmam_fun = bedotbc 
c
	return
	end function wmam_fun