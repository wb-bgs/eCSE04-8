cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine cpt_dat_vals_p
c		Vincent Lesur 09/02/2005
c
c       Parallel interface for cpt_dat_vals.f
c
c       Called: cpt_dat_vals_p
c
c       input:
c         shdeg         max SH degree value
c         nb            Number of base functions
c         nlocpts       number of data+sampling points local to rank
c         nlocdatpts    number of data points local to rank
c         bc            base coefficients
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p(shdeg, nb, nlocpts, nlocdatpts,
     >                            bc, xyzf)
c
        use cudafor, only : cudaDeviceSynchronize
        use kernels
c
        implicit none
c
        integer shdeg, nb
        integer nlocpts, nlocdatpts
c
        real*8 bc(1:nb)
        real*8 xyzf(1:nlocpts)
c
        integer istat
c
c
        call init_cpt_device_arrays(nb, bc)
c
c
c  Offload data points
c          
c
#if defined(CUDA_KERNEL_LOOP)
c
        call cpt_dat_loop(shdeg, 1, nlocdatpts)
c      
#else       
c
        call cpt_dat_kernel <<< get_nblocks_dat(), get_nthreads() >>>
     >    (shdeg, 1, nlocdatpts, 0)
c
#endif
c    
        istat = cudaDeviceSynchronize()
c
#if defined(CUDA_DEBUG)
        call check_for_cuda_error('cpt_dat')
#endif
c
c
c  Offload sampling points
c  Workload is higher due to cpt_sam_kernel_loop()/cpt_sam_kernel() calling XYZsph_bi0_sample()
c
c
#if defined(CUDA_KERNEL_LOOP)
c
        call cpt_sam_loop(shdeg, nlocdatpts+1, nlocpts)
c
#else
c
        call cpt_sam_kernel <<< get_nblocks_sam(), get_nthreads() >>>
     >    (shdeg, nlocdatpts+1, nlocpts, nlocdatpts)
c
#endif
c    
        istat = cudaDeviceSynchronize()
c
#if defined(CUDA_DEBUG)
        call check_for_cuda_error('cpt_sam')
#endif
c
c
        call get_cpt_device_arrays(nlocpts, xyzf)
c
c
        end subroutine cpt_dat_vals_p