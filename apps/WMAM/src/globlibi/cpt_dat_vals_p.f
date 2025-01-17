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
c         nd            Space dimension
c         nlocpts       number of data+sampling points local to rank
c         nlocdatpts    number of data points local to rank
c         d2a           pre-computed array used by mk_lf_dlf()
c         (d)dlf        pre-allocated arrays computed by mk_lf_dlf() and
c                       used within XYZsph_bi0
c         bc            base coefficients
c         ppos(nd+1,*)  point position in nd
c
c       output:
c         XYZF(*)       X,Y,Z or F value at point position
c        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine cpt_dat_vals_p(shdeg, nb, nd, nlocpts, nlocdatpts,
     >                            d2a, dlf, ddlf, bc, ppos, xyzf)
c
        implicit none
c
        integer, unified :: shdeg, nb, nd, nlocpts, nlocdatpts
        real(8), unified :: d2a(0:shdeg)
        real(8), unified :: dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real(8), unified :: bc(1:nb), ppos(nd+1,nlocpts)
        real(8), unified :: xyzf(1:nlocpts)
c
        integer n_threads, n_blocks
c
c
        n_threads = 128
c
        n_blocks = nlocdatpts / n_threads
        call cpt_dat_vals_p_dat<<<n_blocks,n_threads>>>
     >    (shdeg, nb, nd, nlocpts, nlocdatpts,
     >     d2a, dlf, ddlf, bc, ppos, xyzf)
c
        n_blocks = (nlocpts - nlocdatpts) / n_threads
        call cpt_dat_vals_p_smp<<<n_blocks,n_threads>>>
     >    (shdeg, nb, nd, nlocpts, nlocdatpts,
     >     d2a, dlf, ddlf, bc, ppos, xyzf)
c
        end subroutine cpt_dat_vals_p
c
c
c
        attributes(global)
     >  subroutine cpt_dat_vals_p_dat(shdeg, nb, nd,
     >                                nlocpts, nlocdatpts,
     >                                d2a, dlf, ddlf,
     >                                bc, ppos, xyzf)
c
        use XYZsph_bi0
c
        implicit none
c
        integer shdeg, nb, nd
        integer nlocpts, nlocdatpts
        real(8) d2a(0:shdeg)
        real(8) dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real(8) bc(1:nb), ppos(nd+1,nlocpts)
        real(8) xyzf(1:nlocpts)
c
        real(8), parameter :: RAG = 6371.2d0
        real(8), parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
        integer i
c
        real(8) p1, p2, ra
        real(8) bex, bey, bez
c
c
        i = (blockidx%x-1) * blockdim%x
     >      + threadidx%x
c
        if (i .ge. 1 .and.
     >      i .le. nlocdatpts) then
c
          p1 = ppos(1,i)*D2R
          p2 = ppos(2,i)*D2R
          ra = RAG / ppos(3,i)
c
          bex = ppos(5,i)
          bey = ppos(6,i)
          bez = ppos(7,i)
c
          xyzf(i) = XYZsph_bi0_fun(shdeg, nb, d2a,
     >                             dlf, ddlf,
     >                             bc, p1, p2, ra,
     >                             bex, bey, bez)
c
        endif
c
        end subroutine cpt_dat_vals_p_dat
   
   
        attributes(global)
     >  subroutine cpt_dat_vals_p_smp(shdeg, nb, nd,
     >                                nlocpts, nlocdatpts,
     >                                d2a, dlf, ddlf,
     >                                bc, ppos, xyzf)
c
        use XYZsph_bi0
c 
        implicit none
c
        integer shdeg, nb, nd
        integer nlocpts, nlocdatpts
c
        real(8) d2a(0:shdeg)
        real(8) dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real(8) bc(1:nb), ppos(nd+1,nlocpts)
        real(8) xyzf(1:nlocpts)
c
        real(8), parameter :: RAG = 6371.2d0
        real(8), parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
        integer i
c
        real(8) p1, p2, ra
        real(8) bex, bey, bez
c
c   
        i = (blockidx%x-1) * blockdim%x
     >      + threadidx%x
     >      + nlocdatpts
c
        if (i .gt. nlocdatpts .and. 
     >      i .le. nlocpts) then

          p1 = ppos(1,i)*D2R
          p2 = ppos(2,i)*D2R
          ra = RAG / ppos(3,i)
c
          bex = ppos(5,i)
          bey = ppos(6,i)
          bez = ppos(7,i)
c
          call XYZsph_bi0_sample(shdeg, nb, d2a,
     >                           dlf, ddlf, bc,
     >                           p1, p2, ra, 
     >                           bex, bey, bez)
c
          xyzf(i) = XYZsph_bi0_fun(shdeg, nb, d2a,
     >                             dlf, ddlf, bc,
     >                             p1, p2, ra,
     >                             bex, bey, bez)
c
        endif
c
        end subroutine cpt_dat_vals_p_smp