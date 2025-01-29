        module kernels
c
          implicit none
c
c
          public
     >      allocate_device_arrays,
     >      deallocate_device_arrays,
     >      init_nblocks_nthreads,
     >      init_device_arrays,
     >      init_cpt_device_arrays,
     >      init_ssqgh_device_arrays,
     >      get_nblocks_dat,
     >      get_nblocks_sam,
     >      get_nthreads,
     >      get_cpt_device_arrays,
     >      get_ssqgh_device_arrays, 
     >      write_nblocks_nthreads,
     >      cpt_dat_kernel,
     >      cpt_sam_kernel,
     >      ssqgh_dat_kernel,
     >      ssqgh_sam_kernel,
#if defined(CUDA_KERNEL_LOOP)
     >      cpt_dat_loop,
     >      cpt_sam_loop,
     >      cpt_sam_loop_kernel,
     >      ssqgh_dat_loop,
     >      ssqgh_sam_loop,
     >      ssqgh_sam_loop_kernel,
#endif
     >      check_for_cuda_error
c
          private
     >      XYZsph_bi0_sample,
     >      XYZsph_bi0_cpt,
     >      XYZsph_bi0_ssqgh,
     >      mk_lf_dlf
c
c
          private
            real(8), parameter :: RAG = 6371.2d0
            real(8), parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
            integer nblocks_dat, nblocks_sam, nthreads
c
            real(8), allocatable, device :: d_d2a(:)
            real(8), allocatable, device :: d_ppos(:,:)
            real(8), allocatable, device :: d_cov(:)
            integer, allocatable, device :: d_jcov(:)       
c
            real(8), allocatable, device :: d_bc(:)
            real(8), allocatable, device :: d_ddat(:)
            real(8), allocatable, device :: d_xyzf(:) 
            real(8), allocatable, device :: d_gj(:)
            real(8), allocatable, device :: d_dh(:)   
c
c
        contains
c
c
          attributes(host)
     >    subroutine allocate_device_arrays(shdeg, nb, nd, nlocpts)
c
          implicit none
c
          integer, value :: shdeg, nb, nd, nlocpts
c
c
          allocate(d_d2a(0:shdeg))
          allocate(d_ppos(1:nd+1,1:nlocpts))
          allocate(d_cov(1:nlocpts))
          allocate(d_jcov(1:nlocpts+2))
c
          allocate(d_bc(1:nb))
          allocate(d_ddat(1:nlocpts))
          allocate(d_xyzf(1:nlocpts))
          allocate(d_gj(1:nb))
          allocate(d_dh(1:nb))
c
          end subroutine allocate_device_arrays
c
c
c
          attributes(host)
     >    subroutine deallocate_device_arrays()
c
          implicit none
c
c
          deallocate(d_d2a, d_ppos)
          deallocate(d_cov, d_jcov)
c
          deallocate(d_bc, d_ddat, d_xyzf)
          deallocate(d_gj, d_dh)
c     
          end subroutine deallocate_device_arrays
c
c
c
          attributes(host)
     >    subroutine init_nblocks_nthreads(cuda_nblocks_dat,
     >                                     cuda_nblocks_sam,
     >                                     cuda_nthreads,
     >                                     nlocdatpts, nlocsampts)
c
          use cudafor
c
          implicit none
c
          integer cuda_nblocks_dat, cuda_nblocks_sam
          integer cuda_nthreads
          integer nlocdatpts, nlocsampts
c         
          type(cudaDeviceProp) :: cuda_prop
c
          integer ierr, maxblocks, maxthreads
c
c
          ierr = cudaGetDeviceProperties(cuda_prop, 0)
c
          maxblocks = cuda_prop%multiProcessorCount
     >              * cuda_prop%maxBlocksPerMultiProcessor
c
          maxthreads = cuda_prop%maxThreadsPerBlock
c
          if (cuda_nthreads .le. 0 .or.
     >        cuda_nthreads .gt. maxthreads) then
            nthreads = 128
          else
            nthreads = cuda_nthreads
          endif
c
          if (cuda_nblocks_dat .le. 0 .or.
     >        cuda_nblocks_dat .gt. maxblocks) then
            nblocks_dat = nlocdatpts / nthreads
            if (MOD(nlocdatpts, nthreads) .gt. 0) then
              nblocks_dat = nblocks_dat + 1
            endif
          else
            nblocks_dat = cuda_nblocks_dat
          endif
c
          if (cuda_nblocks_sam .le. 0 .or.
     >        cuda_nblocks_sam .gt. maxblocks) then
            nblocks_sam = nlocsampts / nthreads
            if (MOD(nlocsampts, nthreads) .gt. 0) then
              nblocks_sam = nblocks_sam + 1
            endif
          else
            nblocks_sam = cuda_nblocks_sam
          endif
c
          end subroutine init_nblocks_nthreads
c
c
c
          attributes(host)
     >    subroutine init_device_arrays(shdeg, nd, nlocpts,
     >                                  d2a, ppos, cov, jcov)
c          
          implicit none
c
          integer shdeg, nd
          integer nlocpts, nlocdatpts
c
          real*8 d2a(0:shdeg)
          real*8 ppos(1:nd+1,1:nlocpts)
          real*8 cov(1:nlocpts)
          integer jcov(1:nlocpts+2)
c
c
          d_d2a  = d2a
          d_ppos = ppos
          d_cov  = cov
          d_jcov = jcov
c
          end subroutine init_device_arrays
c
c
c
          attributes(host)
     >    subroutine init_cpt_device_arrays(nb, bc)
c          
          implicit none
c
          integer, value :: nb
c
          real*8 bc(1:nb)
c
c
          d_bc = bc
c
          end subroutine init_cpt_device_arrays
c
c
c
          attributes(host)
     >    subroutine init_ssqgh_device_arrays(nb, nlocpts,
     >                                        bc, ddat, xyzf)
c          
          implicit none
c
          integer, value :: nb, nlocpts
c
          real*8 bc(1:nb)
          real*8 ddat(1:nlocpts)
          real*8 xyzf(1:nlocpts)
c
c
          d_bc   = bc
          d_ddat = ddat
          d_xyzf = xyzf
c
          d_gj   = 0.0d0
          d_dh   = 0.0d0
c
          end subroutine init_ssqgh_device_arrays
c
c
c
          integer function get_nblocks_dat()
c
          implicit none
c
c
          get_nblocks_dat = nblocks_dat
c
          end function get_nblocks_dat
c
c
c
          integer function get_nblocks_sam()
c
          implicit none
c
c
          get_nblocks_sam = nblocks_sam
c
          end function get_nblocks_sam
c
c
c
          integer function get_nthreads()
c
          implicit none
c
c
          get_nthreads = nthreads
c
          end function get_nthreads
c
c
c
          attributes(host)
     >    subroutine get_cpt_device_arrays(nlocpts, xyzf)
c          
          implicit none
c
          integer, value :: nlocpts
c
          real*8 xyzf(1:nlocpts)
c
c
          xyzf = d_xyzf
c
          end subroutine get_cpt_device_arrays
c
c           
c
          attributes(host)              
     >    subroutine get_ssqgh_device_arrays(nb, gj, dh)
c          
          implicit none
c
          integer, value :: nb
c
          real*8 gj(1:nb)
          real*8 dh(1:nb)
c
c
          gj = d_gj
          dh = d_dh
c
          end subroutine get_ssqgh_device_arrays
c
c
c
          attributes(host)
     >    subroutine write_nblocks_nthreads()
c
          implicit none
c
          write(*,*) ''
#if defined(CUDA_KERNEL_LOOP) && !defined(CUDA_KERNEL_LOOP_USER)
          write(*,*) 'cuda_nblocks_dat: *'
          write(*,*) 'cuda_nblocks_sam: *'
          write(*,*) 'cuda_nthreads: *'
#else
          write(*,*) 'cuda_nblocks_dat: ', nblocks_dat
          write(*,*) 'cuda_nblocks_sam: ', nblocks_sam
          write(*,*) 'cuda_nthreads: ', nthreads
#endif
          write(*,*) ''
          write(*,*) ''
c
          end subroutine write_nblocks_nthreads
c
c
c
          attributes(global)
     >    subroutine cpt_dat_kernel(shdeg, imin, imax, ioffset)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg
          integer, value :: imin, imax, ioffset
c
          integer i
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
c
c
          i = (blockidx%x-1) * blockdim%x
     >      + threadidx%x
     >      + ioffset
c
          if (i .ge. imin .and.
     >        i .le. imax) then
c
            p1 = d_ppos(1,i)*D2R
            p2 = d_ppos(2,i)*D2R
            ra = RAG / d_ppos(3,i)
c
            bex = d_ppos(5,i)
            bey = d_ppos(6,i)
            bez = d_ppos(7,i)
c
            d_xyzf(i) = XYZsph_bi0_cpt(shdeg, p1, p2, ra,
     >                                 bex, bey, bez)
c
          endif
c
          end subroutine cpt_dat_kernel
c
c
c
          attributes(global)
     >    subroutine cpt_sam_kernel(shdeg, imin, imax, ioffset)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg
          integer, value :: imin, imax, ioffset
c
          integer i
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
c
c
          i = (blockidx%x-1) * blockdim%x
     >      + threadidx%x
     >      + ioffset
c
          if (i .ge. imin .and.
     >        i .le. imax) then
c
            p1 = d_ppos(1,i)*D2R
            p2 = d_ppos(2,i)*D2R
            ra = RAG / d_ppos(3,i)
c
            bex = d_ppos(5,i)
            bey = d_ppos(6,i)
            bez = d_ppos(7,i)
c
            call XYZsph_bi0_sample(shdeg, p1, p2, ra,
     >                             bex, bey, bez)
c
            d_xyzf(i) = XYZsph_bi0_cpt(shdeg, p1, p2, ra,
     >                                 bex, bey, bez)
c
          endif
c
          end subroutine cpt_sam_kernel
c
c
c
          attributes(global)
     >    subroutine ssqgh_dat_kernel(shdeg, imin, imax, ioffset)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg
          integer, value :: imin, imax, ioffset
c
          integer i
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
          real(8) dw_gj, dw_dh
c
c
          i = (blockidx%x-1) * blockdim%x
     >      + threadidx%x
     >      + ioffset
c
          if (i .ge. imin .and. i .le. imax) then
c      
            p1 = d_ppos(1,i)*D2R
            p2 = d_ppos(2,i)*D2R
            ra = RAG / d_ppos(3,i)
c
            bex = d_ppos(5,i)
            bey = d_ppos(6,i)
            bez = d_ppos(7,i)
c
            dw_dh = 2.d0*(1.d0/d_cov(d_jcov(i)))
            dw_gj = dw_dh*(d_ddat(i)-d_xyzf(i))
c
            call XYZsph_bi0_ssqgh(shdeg, p1, p2, ra,
     >                            bex, bey, bez,
     >                            dw_gj, dw_dh)
c
          endif
c
          end subroutine ssqgh_dat_kernel
c
c
c
          attributes(global)
     >    subroutine ssqgh_sam_kernel(shdeg, imin, imax, ioffset)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg
          integer, value :: imin, imax, ioffset
c
          integer i
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
          real(8) dw_gj, dw_dh
c
c
          i = (blockidx%x-1) * blockdim%x
     >      + threadidx%x
     >      + ioffset
c
          if (i .ge. imin .and. i .le. imax) then
c      
            p1 = d_ppos(1,i)*D2R
            p2 = d_ppos(2,i)*D2R
            ra = RAG / d_ppos(3,i)
c
            bex = d_ppos(5,i)
            bey = d_ppos(6,i)
            bez = d_ppos(7,i)
c
            dw_dh = 2.d0*(1.d0/d_cov(d_jcov(i)))
            dw_gj = dw_dh*(d_ddat(i)-d_xyzf(i))
c
            call XYZsph_bi0_sample(shdeg, p1, p2, ra,
     >                             bex, bey, bez)
c
            call XYZsph_bi0_ssqgh(shdeg, p1, p2, ra,
     >                            bex, bey, bez,
     >                            dw_gj, dw_dh)
c
          endif
c
          end subroutine ssqgh_sam_kernel
c
c
c
#if defined(CUDA_KERNEL_LOOP)
c
c
          attributes(host)
     >    subroutine cpt_dat_loop(shdeg, imin, imax)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg, imin, imax
c
          integer i
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
c
c
#if defined(CUDA_KERNEL_LOOP_USER)
!$cuf kernel do <<< get_nblocks_dat(), get_nthreads() >>>
#else
!$cuf kernel do <<< *, * >>>
#endif
          do i = imin,imax
c
            p1 = d_ppos(1,i)*D2R
            p2 = d_ppos(2,i)*D2R
            ra = RAG / d_ppos(3,i)
c
            bex = d_ppos(5,i)
            bey = d_ppos(6,i)
            bez = d_ppos(7,i)
c
            d_xyzf(i) = XYZsph_bi0_cpt(shdeg, p1, p2, ra,
     >                                 bex, bey, bez)
c
          enddo
c
          end subroutine cpt_dat_loop
c
c
c
          attributes(host)
     >    subroutine cpt_sam_loop(shdeg, imin, imax)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg, imin, imax
c
          integer i
c
c
#if defined(CUDA_KERNEL_LOOP_USER)
!$cuf kernel do <<< get_nblocks_sam(), get_nthreads() >>>
#else
!$cuf kernel do <<< *, * >>>
#endif
          do i = imin,imax
c
            call cpt_sam_loop_kernel(shdeg, i)
c
          enddo
c
          end subroutine cpt_sam_loop
c
c
c
          attributes(device)
     >    subroutine cpt_sam_loop_kernel(shdeg, ip)
c
          implicit none
c
          integer, value :: shdeg, ip
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
c
          p1 = d_ppos(1,ip)*D2R
          p2 = d_ppos(2,ip)*D2R
          ra = RAG / d_ppos(3,ip)
c
          bex = d_ppos(5,ip)
          bey = d_ppos(6,ip)
          bez = d_ppos(7,ip)
c
          call XYZsph_bi0_sample(shdeg, p1, p2, ra,
     >                           bex, bey, bez)
c
          d_xyzf(ip) = XYZsph_bi0_cpt(shdeg, p1, p2, ra,
     >                                bex, bey, bez)
c
          end subroutine cpt_sam_loop_kernel
c
c
c
          attributes(host)
     >    subroutine ssqgh_dat_loop(shdeg, imin, imax)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg, imin, imax
c
          integer i
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
          real(8) dw_gj, dw_dh
c
c
#if defined(CUDA_KERNEL_LOOP_USER)
!$cuf kernel do <<< get_nblocks_dat(), get_nthreads() >>>
#else
!$cuf kernel do <<< *, * >>>
#endif
          do i = imin,imax
c
            p1 = d_ppos(1,i)*D2R
            p2 = d_ppos(2,i)*D2R
            ra = RAG / d_ppos(3,i)
c
            bex = d_ppos(5,i)
            bey = d_ppos(6,i)
            bez = d_ppos(7,i)
c
            dw_dh = 2.d0*(1.d0/d_cov(d_jcov(i)))
            dw_gj = dw_dh*(d_ddat(i)-d_xyzf(i))
c
            call XYZsph_bi0_ssqgh(shdeg, p1, p2, ra,
     >                            bex, bey, bez,
     >                            dw_gj, dw_dh)
c
          enddo
c
          end subroutine ssqgh_dat_loop
c
c
c
          attributes(host)
     >    subroutine ssqgh_sam_loop(shdeg, imin, imax)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg, imin, imax
c
          integer i
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
          real(8) dw_gj, dw_dh
c
c
#if defined(CUDA_KERNEL_LOOP_USER)
!$cuf kernel do <<< get_nblocks_sam(), get_nthreads() >>>
#else
!$cuf kernel do <<< *, * >>>
#endif
          do i = imin,imax
c
            call cpt_ssqgh_loop_kernel(shdeg, i)
c
          enddo
c
          end subroutine ssqgh_sam_loop
c
c
c
          attributes(device)
     >    subroutine cpt_ssqgh_loop_kernel(shdeg, ip)
c
          implicit none
c
          integer, value :: shdeg, ip
c
          real(8) p1, p2, ra
          real(8) bex, bey, bez
          real(8) dw_gj, dw_dh
c
c
          p1 = d_ppos(1,ip)*D2R
          p2 = d_ppos(2,ip)*D2R
          ra = RAG / d_ppos(3,ip)
c
          bex = d_ppos(5,ip)
          bey = d_ppos(6,ip)
          bez = d_ppos(7,ip)
c
          dw_dh = 2.d0*(1.d0/d_cov(d_jcov(ip)))
          dw_gj = dw_dh*(d_ddat(ip)-d_xyzf(ip))
c
          call XYZsph_bi0_sample(shdeg, p1, p2, ra,
     >                           bex, bey, bez)
c
          call XYZsph_bi0_ssqgh(shdeg, p1, p2, ra,
     >                          bex, bey, bez,
     >                          dw_gj, dw_dh)
c
          end subroutine cpt_ssqgh_loop_kernel
c
c
c  end of <#if defined(CUDA_KERNEL_LOOP)> clause
#endif
c
c
c
          attributes(host)
     >    subroutine check_for_cuda_error(kernel_name)
c
          use cudafor
c
          implicit none
c
          include 'mpif.h'
c
          character(*) kernel_name
c
          integer rank, ierr
c
c
          call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
          ierr = cudaGetLastError()
          if (ierr .gt. 0) then
            write(*,*) rank,
     >        ': Error, ', kernel_name, ' kernel failure: ',
     >        ierr, ', ', cudaGetErrorString(ierr)
          endif
c
          end subroutine check_for_cuda_error
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine XYZsph_bi0_sample
c		
c       Compute dx = bx.bc, dy = by.bc and dz = bz.bc.
c       This subroutine is called for sample points only.
c             
c       input:
c          shdeg  INTEGER       max SH degree value
c          p1     REAL(8)       co-latitude
c          p2     REAL(8)       longitude
c          ra     REAL(8)       radius
c
c       input/output:
c          bex    REAL(8)       x component of magnetic field
c          bey    REAL(8)       y component of magnetic field
c          bez    REAL(8)       z component of magnetic field
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          attributes(device)
     >    subroutine XYZsph_bi0_sample(shdeg, p1, p2, ra,
     >                                 bex, bey, bez)
c
          implicit none
c
          integer, value :: shdeg
c
          real(8), value :: p1, p2, ra
          real(8) bex, bey, bez
c
          integer nu, il, im, ik
          real(8) rc, rs
          real(8) ds, dc, dr, dw
          real(8) bx, by, bz
          real(8) bxp1, byp1, bzp1
          real(8) bc_nu, bc_nup1
c
          real(8) dx, dy, dz, dd
          real(8) dxbey, dxbez
          real(8) dybex, dybez
          real(8) dzbex, dzbey
c
          real(8) xy_c, xz_c
          real(8) yx_c, yz_c
          real(8) zx_c, zy_c
c
          real(8) bex2, bey2, bez2
c
          real(8) dlf(1:shdeg+1)
          real(8) ddlf(1:shdeg+1)
c        
c 
          dx = 0.0d0
          dy = 0.0d0 
          dz = 0.0d0 
c
          rc = dcos(p1)
          rs = dsin(p1)
c
          call mk_lf_dlf(0, shdeg, rs, rc, dlf, ddlf)
          do il = 1,shdeg
            dr = ra**(il+2)
c          
            ik = il+1
c          
            bx =  ddlf(ik) * dr
            bz = -dlf(ik)  * dr * dble(ik)
c
            dx = dx + bx * d_bc(il)
            dz = dz + bz * d_bc(il)
          enddo
c
c
          nu = shdeg + 1
          do im = 1,shdeg
            call mk_lf_dlf(im, shdeg, rs, rc, dlf, ddlf)
            dc = dcos(im*p2)
            ds = dsin(im*p2)
            do il = im,shdeg
              bc_nu   = d_bc(nu)
              bc_nup1 = d_bc(nu+1)

              ik = il-im+1
              dr = ra**(il+2)

              dw   = ddlf(ik) * dr
              bx   = dc * bc_nu
              bxp1 = ds * bc_nup1
              dx   = dx + dw*(bx + bxp1)

              dw   = (dlf(ik)/rs) * dr * dble(im)
              by   =  ds * bc_nu
              byp1 = -dc * bc_nup1
              dy   =  dy + dw*(by + byp1)
            
              dw   =  dlf(ik) * dr * dble(il+1)
              bz   = -dc * bc_nu
              bzp1 = -ds * bc_nup1
              dz   =  dz + dw*(bz + bzp1)

              nu = nu + 2
            enddo
          enddo
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
c	
          bex = bex2
          bey = bey2
          bez = bez2
c
          end subroutine XYZsph_bi0_sample
c     
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	function XYZsph_bi0_cpt
c		
c       Computes the dot product of 'be' and 'bc', avoiding the need
c       to store the entire contents of the 'be' coefficient array.
c             
c       input:
c          shdeg  INTEGER     max SH degree value
c          p1     REAL(8)     co-latitude
c          p2     REAL(8)     longitude
c          ra     REAL(8)     radius
c          bex    REAL(8)     x component of magnetic field
c          bey    REAL(8)     y component of magnetic field
c          bez    REAL(8)     z component of magnetic field
c
c       output:
c          XYZsph_bi0_cpt  REAL(8)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          attributes(device)
     >    real(8)
     >    function XYZsph_bi0_cpt(shdeg, p1, p2, ra,
     >                            bex, bey, bez)
c
          implicit none
c
          integer, value :: shdeg
c
          real(8), value :: p1, p2, ra
          real(8), value :: bex, bey, bez
c
          integer nu, il, im, ik
          real(8) rc, rs
          real(8) ds, dc, dr, dw
          real(8) bx, by, bz
          real(8) bxp1, byp1, bzp1
c
          real(8) dlf(1:shdeg+1)
          real(8) ddlf(1:shdeg+1)
c 
c
          XYZsph_bi0_cpt = 0.0d0 
c        
          rc = dcos(p1)
          rs = dsin(p1)
c
          call mk_lf_dlf(0, shdeg, rs, rc, dlf, ddlf)
          do il = 1,shdeg
            dr = ra**(il+2)
c          
            ik = il+1
c          
            bx =  ddlf(ik) * dr
            bz = -dlf(ik)  * dr * dble(ik)
c
            XYZsph_bi0_cpt = XYZsph_bi0_cpt
     >                     + (bex*bx + bez*bz) * d_bc(il)
          enddo
c
c
          nu = shdeg + 1
          do im = 1,shdeg
            call mk_lf_dlf(im, shdeg, rs, rc, dlf, ddlf)
            dc = dcos(im*p2)
            ds = dsin(im*p2)
            do il = im,shdeg
              ik = il-im+1
              dr = ra**(il+2)
c
              dw   = ddlf(ik) * dr * bex
              bx   = dw*dc
              bxp1 = dw*ds
c            
              dw   = (dlf(ik)/rs) * dr * dble(im) * bey
              by   =  dw*ds
              byp1 = -dw*dc
c            
              dw   =  dlf(ik) * dr * dble(il+1) * bez
              bz   = -dw*dc
              bzp1 = -dw*ds
c
              XYZsph_bi0_cpt = XYZsph_bi0_cpt
     >                       + (bx + by + bz) * d_bc(nu)
     >                       + (bxp1 + byp1 + bzp1) * d_bc(nu+1)
c
              nu = nu + 2
            enddo
          enddo
c
          end function XYZsph_bi0_cpt
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine XYZsph_bi0_ssqgh
c		
c       Computes the gradient of the weighted sum of squares (gj) and
c       the diagonal of the Hessian (hj). 
c
c       input:
c          shdeg    INTEGER     max SH degree value
c          p1       REAL(8)     co-latitude
c          p2       REAL(8)     longitude
c          ra       REAL(8)     radius
c          bex      REAL(8)     x component of magnetic field
c          bey      REAL(8)     y component of magnetic field
c          bez      REAL(8)     z component of magnetic field 
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          attributes(device)
     >    subroutine XYZsph_bi0_ssqgh(shdeg, p1, p2, ra,
     >                                bex, bey, bez,
     >                                dw_gj, dw_dh)
c
          use cudafor
c
          implicit none
c
          integer, value :: shdeg      
c
          real(8), value :: p1, p2, ra
          real(8), value :: bex, bey, bez 
          real(8), value :: dw_gj, dw_dh
c
          integer nu, il, im, ik, istat
          real(8) rc, rs
          real(8) ds, dc, dr, dw
          real(8) bx, by, bz
          real(8) bxp1, byp1, bzp1
          real(8) be, bep1
c
          real(8) dlf(1:shdeg+1)
          real(8) ddlf(1:shdeg+1)
c 
c
          rc = dcos(p1)
          rs = dsin(p1)
c
          call mk_lf_dlf(0, shdeg, rs, rc, dlf, ddlf)
          do il = 1,shdeg
            dr = ra**(il+2)
c          
            ik = il+1
c          
            bx =  ddlf(ik) * dr
            bz = -dlf(ik)  * dr * dble(ik)
c
            be = (bex*bx + bez*bz)
c
            istat = atomicadd(d_gj(il), dw_gj*be)
            istat = atomicadd(d_dh(il), dw_dh*(be**2))
          enddo
c
c
          nu = shdeg + 1
          do im = 1,shdeg
            call mk_lf_dlf(im, shdeg, rs, rc, dlf, ddlf)
            dc = dcos(im*p2)
            ds = dsin(im*p2)
            do il = im,shdeg
              ik = il-im+1
              dr = ra**(il+2)
c
              dw   = ddlf(ik) * dr
              bx   = dw*dc
              bxp1 = dw*ds
c            
              dw   = (dlf(ik)/rs) * dr * dble(im)
              by   =  dw*ds
              byp1 = -dw*dc
c            
              dw   =  dlf(ik) * dr * dble(il+1)
              bz   = -dw*dc
              bzp1 = -dw*ds
c
              be   = bex*bx + bey*by + bez*bz
              bep1 = bex*bxp1 + bey*byp1 + bez*bzp1
c
              istat = atomicadd(d_gj(nu), dw_gj*be)
              istat = atomicadd(d_gj(nu+1), dw_gj*bep1)
c
              istat = atomicadd(d_dh(nu), dw_dh*(be**2))
              istat = atomicadd(d_dh(nu+1), dw_dh*(bep1**2))
c
              nu = nu + 2
            enddo
          enddo
c
          end subroutine XYZsph_bi0_ssqgh
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
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
c         im/shdeg      order and degree max
c         rc/rs         cos(colatitude)/sin(colatitude)
c       output:
c         dlf           legendre function from im to nl
c         ddlf          derivative of legendre function from im to nl
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          attributes(device)
     >    subroutine mk_lf_dlf(im, shdeg, rs, rc, dlf, ddlf)
c
          implicit none
c
          integer im, shdeg
c
          real(8) rs, rc
c
          real(8) dlf(1:shdeg+1)
          real(8) ddlf(1:shdeg+1)
c
          integer d0, il, jl
          real(8) d1, d2
          real(8) dalpha, dbeta
c
c
c  Initialise dlf array
          d1 = rs
          if (d1 .ne. 0.0d0) then
            d1 = d1**im
          else
            if (im .eq. 0) d1 = 1.d0
          endif
c
          dlf(1) = d1*d_d2a(im)
          dlf(2) = dlf(1) * rc * dsqrt(dble(2*im+1))
c
          do il = 2,shdeg-im
            d0 = il+2*im
            d1 = dble((il-1) * (d0-1))
            d2 = dble(il * d0)
            dbeta = dsqrt(d1/d2)*dlf(il-1)
c
            d1 = dble(2*(il+im)-1)
            dalpha = (d1/dsqrt(d2))*dlf(il)*rc
c
            dlf(il+1) = dalpha - dbeta
          enddo
c
c
c  Initialise ddlf array
          if (rs .eq. 0.0d0) then
c
            if (im .ne. 1) then 
              ddlf(1:shdeg+1) = 0.0d0
            else
              ddlf(1) = -1.0d0
              ddlf(2) = -dsqrt(3.0d0) 
              do il = 3,shdeg+1
                jl = (il**2)-1
                d1 = (2*il-1)/dsqrt(dble(jl))
                d2 = dsqrt(dble((il-1)**2-1)/jl)
                ddlf(il) = d1*ddlf(il-1) - d2*ddlf(il-2)
              enddo
            endif
c
          else
c
            do il = shdeg,im+1,-1
              jl = il - im + 1
              d1 = dsqrt(dble((il-im)*(il+im)))
              d2 = dble(il) 
              ddlf(jl) = (d2*rc*dlf(jl) - d1*dlf(jl-1)) / rs
            enddo
            ddlf(1) = im * rc * dlf(1)/rs
c
          endif
c
          end subroutine mk_lf_dlf
c
c
        end module kernels