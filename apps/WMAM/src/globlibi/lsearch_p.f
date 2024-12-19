ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c          subroutine LSEARCH_P
c                           vincent Lesur 14/06/2006
c
c      Modified 18.08.2011 to call cptstd_d in place of cptstd as
c      none diagonal covariance matrix element are not expected
c      for Large system (V.Lesur)
c
c      2017-09-14 (bham) Replaced expression (now) on line 212:
c        (stp-st(im1))/st(im1)
c       with:
c        (stp-st(im1))/max(abs(st(im1)),abs(stp),epss)
c       because were having problems with `st(im1)` being zero.
c      15.07.2011 (v.Lesur) update input list of cptstd
c      14.07.2011 FM & FDAMP removed from parameter list
c      10.3.2011 sign of ds changed (V.Lesur)
c
c     Look for the minimun of a non-linear function in the direction DS
c     Parallel version
c
c     input:
c         iunit         unit to write outputs
c         itm           Maximum number of iterations
c         shdeg         max SH degree value
c         nb            Number or base function to use
c         nd            space dimension
c         npts          Total number of points (data + sampling) for all ranks
c         nlocpts       Total number of points for this rank
c         nlocdatpts    number of data points assigned to rank
c         d2a           pre-computed array used by mk_lf_dlf()
c         dra           pre-allocated array used within XYZsph_bi0
c         dlf           "
c         ddlf          "
c         bc            Estimate of Base function coefficients
c         ppos          data point position in ndD
c         ddat          data values
c         src_stat      MPI gradient search status
c         dl(3)         control lsearch process & damping
c         cov           covariance matrix in SLAP Column format
c         jcov          Integer vector describing cov format
c         std           STD value for given BC
c         ds            gradient direction
c
c       output:
c         stp           recommended step in direction ds(*)
c         std           STD value for given BC+stp*DS
c         xyzf          Forward modelling for given BC+stp*DS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine lsearch_p(iunit, itm, shdeg, nb, nd,
     >                       npts, nlocpts, nlocdatpts,
     >                       d2a, dra, dlf, ddlf,
     >                       bc, ppos, ddat,
     >                       src_stat, MPI_SEARCH_STATUS,
     >                       dl, cov, jcov,
     >                       std, ds, stp, xyzf)
c
        implicit none
c
        include 'mpif.h'
        include 'mpi_status_types.h'
c
        integer iunit, itm, shdeg, nb, nd
        integer npts, nlocpts, nlocdatpts
        real*8 d2a(0:shdeg), dra(1:shdeg)
        real*8 dlf(1:shdeg+1), ddlf(1:shdeg+1)
        real*8 bc(1:nb)
        real*8 ppos(1:nd+1,1:nlocpts)
        real*8 ddat(1:nlocpts)
        type(search_status) src_stat
        integer MPI_SEARCH_STATUS
        real*8 dl(1:3), cov(1:nlocpts)
        integer jcov(1:nlocpts+2)
        real*8 std, ds(1:nb), stp, xyzf(1:nlocpts)
c
        integer i, im1, im2, it
        integer ierr, rank
        real*8 dj(1:3), st(1:3)
        real*8 fct, fctt, fctl, fcts, dd, rt
        real*8 epss, stp1, stp2, gr, gl, numer, denom
        real*8, allocatable :: bcn(:)
        character yon_rf
c
c
c  All defining parallel enviroment
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
 
        src_stat%stp=stp
c
        allocate(bcn(1:nb))
c
c MP: does the all lot
        if (rank .eq. 0) then
c
c MP0-Initial
            src_stat%yon_ct = 'y'
            epss = dl(1)
            fctl = 5.d0
            fcts = 1.d0/fctl
            rt = 4.d0
c
            fct = 2.d0
            if (src_stat%stp .eq. 0.0d0) src_stat%stp = 1.d0
            src_stat%stp = src_stat%stp/fct
c
c MP1-Find the large steps
c
c MP 11- initial (given in input: �std� no need to compute)
            dj(1) = std
            st(1) = 0.0d0
            dj(2) = std
            st(2) = 0.0d0
c
c MP 12- find st(3) such that dj(3) significantly different from dj(1)
            dd = 0.d0
            it = 1
            do while (dd .le. epss .and. it .le. itm)
                src_stat%stp = src_stat%stp*fct
c               write(*,*)'lsearch_p: 1'
                call MPI_BCAST(src_stat, 1, MPI_SEARCH_STATUS,
     >                         0, MPI_COMM_WORLD, ierr)
c
                bcn(1:nb) = bc(1:nb)+src_stat%stp*ds(1:nb)
c
                call cpt_dat_vals_p(shdeg, nb, nd,
     >                              nlocpts, nlocdatpts,
     >                              d2a, dra, dlf, ddlf,
     >                              bcn, ppos, xyzf)
c
                call cptstd_dp(npts, nlocpts,
     >                         cov, jcov, ddat,
     >                         xyzf, std)
c
                dj(3) = std
                dd = dabs(dj(3)-dj(1))/dj(1)
                it = it+1
            enddo

            if (it .gt. itm) then
                write(iunit,*) 'Lsearch: cannot find minimum step'
                src_stat%yon_ct = 'n'
                im1 = 1
            endif
            st(3) = src_stat%stp
c
c MP 13- Find st(i),st(im1),st(im2)such that dj(im1)<dj(i) & dj(im1)<dj(im2) 
            if (src_stat%yon_ct .eq. 'y') then
                it = 1
                if (dj(3) .lt. dj(1)) then
                    im1 = 2
                    i = 3
                    fctt = fct
                    do while (dj(i) .lt. dj(im1) .and. it .le. itm)
                        im1 = i
                        i = mod(im1,3)+1
                        src_stat%stp = src_stat%stp*fctt
                        fctt = fctt*fct
c                       write(*,*)'lsearch_p: 2'
                        call MPI_BCAST(src_stat, 1, MPI_SEARCH_STATUS,
     >                                 0, MPI_COMM_WORLD, ierr)
c
                        bcn(1:nb) = bc(1:nb)+src_stat%stp*ds(1:nb)
c
                        call cpt_dat_vals_p(shdeg, nb, nd,
     >                                      nlocpts, nlocdatpts,
     >                                      d2a, dra, dlf, ddlf,
     >                                      bcn, ppos, xyzf)
c
                        call cptstd_dp(npts, nlocpts,
     >                                 cov, jcov, ddat,
     >                                 xyzf, std)
c
                        dj(i) = std
                        st(i) = src_stat%stp
                        it = it+1
                    enddo

                    im2 = 6-(i+im1)
                    if (it .gt. itm) then
                        write(iunit,*) 
     >                      'Lsearch: cannot increase functional value'
                        src_stat%yon_ct = 'n'
                    endif
c end of <if (dj(3).lt.dj(1))> clause
                else
                    im1 = 3
                    im2 = 1
                    fctt = fct
                    src_stat%stp = -src_stat%stp/fctt
                    do while (dj(im2) .lt. dj(im1) .and. it .le. itm)
                        im1 = im2
                        im2 = mod(im1,3)+1
                        src_stat%stp = src_stat%stp*fctt
                        fctt = fctt*fct
c                       write(*,*)'lsearch_p: 3'
                        call MPI_BCAST(src_stat, 1, MPI_SEARCH_STATUS,
     >                                 0, MPI_COMM_WORLD, ierr)
c
                        bcn(1:nb) = bc(1:nb)+src_stat%stp*ds(1:nb)
c
                        call cpt_dat_vals_p(shdeg, nb, nd,
     >                                      nlocpts, nlocdatpts,
     >                                      d2a, dra, dlf, ddlf,
     >                                      bcn, ppos, xyzf)
c
                        call cptstd_dp(npts, nlocpts,
     >                                 cov, jcov, ddat,
     >                                 xyzf, std)
c
                        dj(im2) = std
                        st(im2) = src_stat%stp
                        it = it+1
                    enddo

                    i = 6-(im2+im1)
                    if (it .gt. itm) then
                        write(iunit,*)
     >                      'Lsearch: cannot increase functional value'
                        src_stat%yon_ct = 'n'
                    endif

c end of <if (dj(3) .ge. dj(1))> clause
                endif

                write(iunit,*)'lsearch start',it
                write(iunit,'(A,3e15.7)')
     >              'lsearch start',dj(i),dj(im1),dj(im2)
                write(iunit,'(A,3e15.7)')
     >              'lsearch start',st(i),st(im1),st(im2)
        
c end of <if (src_stat%yon_ct.eq.'y')> clause
            endif
c
c MP2-Find the zero gradient
            it = 1
            do while (src_stat%yon_ct .eq. 'y')
                it = it+1
                stp1 = st(i)-st(im1)  
                stp2 = st(im1)-st(im2)
c
c MP 20- check for step size
                yon_rf = 'y'
                if(dabs(st(im1)) .gt. epss)then
                    dd = dabs((stp1+stp2)/st(im1))
                else
                    dd = dabs(stp1+stp2)
                endif
                if (dd.lt.epss) yon_rf = 'n'
c
c MP 21- If step size large enough: Find the next step try
                if (yon_rf .eq. 'y') then
c MP  211- check for step interval relative sizes and set default values if needed 
                    fct = dabs(stp1/stp2)
                    if (fct .lt. fcts) then
                        src_stat%stp = st(im1)-(st(im1)-st(im2))/rt
                    elseif (fct .gt. fctl) then
                        src_stat%stp = st(im1)+(st(i)-st(im1))/rt
                    else
c MP  212- steps are similar=> computes gradients (gr*gl<0)
c        stp given by False position algorithm
                        gr = (dj(i)-dj(im1))/stp1
                        gl = (dj(im1)-dj(im2))/stp2
c
                        src_stat%stp = (st(im2)*gr-st(i)*gl)/(gr-gl)
                        src_stat%stp = 0.5d0*(st(im1)+src_stat%stp)
c
c MP  213-check step size and set default values if needed
                        numer = src_stat%stp-st(im1)
                        denom = max(abs(st(im1)),abs(src_stat%stp),epss)
                        if (dabs(numer/denom) .lt. epss) then
                            if (fct .lt. 1.d0) then
                                src_stat%stp = st(im1)
     >                              -(st(im1)-st(im2))/rt
                            else
                                src_stat%stp = st(im1)
     >                              +(st(i)-st(im1))/rt
                            endif
                        endif
                    endif
c
c MP 22-  Find functional values and new im1,im2 & i
c                   write(*,*)'lsearch_p: 4'
                    call MPI_BCAST(src_stat, 1, MPI_SEARCH_STATUS,
     >                             0, MPI_COMM_WORLD, ierr)
c
                    bcn(1:nb) = bc(1:nb)+src_stat%stp*ds(1:nb)
c
                    call cpt_dat_vals_p(shdeg, nb, nd,
     >                                  nlocpts, nlocdatpts,
     >                                  d2a, dra, dlf, ddlf,
     >                                  bcn, ppos, xyzf)
c
                    call cptstd_dp(npts, nlocpts,
     >                             cov, jcov, ddat,
     >                             xyzf, std)
c
                    if (src_stat%stp .gt. st(im1)) then
                        if (std .ge. dj(im1)) then
                            dj(i) = std
                            st(i) = src_stat%stp
                        else
                            dj(im2) = std
                            st(im2) = src_stat%stp
                            im1 = im2
                            im2 = 6-(i+im1)
                        endif
                    else
                        if (std .ge. dj(im1)) then
                            dj(im2) = std
                            st(im2) = src_stat%stp
                        else
                            dj(i) = std
                            st(i) = src_stat%stp
                            im1 = i
                            i = 6-(im2+im1)
                        endif
                    endif
c
c MP 23- Check if maximum iteration 
                    if (it.ge.itm) then
                        src_stat%yon_ct='n'
                        write(iunit,*)'Lsearch: Maximum iteration',it
                        write(iunit,'(3e15.7)')dj(i),dj(im1),dj(im2)
                        write(iunit,'(3e15.7)')st(i),st(im1),st(im2)
                    endif
c
c MP 24- Check if significant improvement still possible
                    dd = max(dj(im2),dj(i))
                    dd = dabs((dd-dj(im1))/dj(im1))
                    if (dd .lt. epss) then 
                        src_stat%yon_ct = 'n'
                        write(iunit,*)
     >                      'Lsearch: No hope of improvement',it
                        write(iunit,'(3e15.7)')dj(i),dj(im1),dj(im2)
                        write(iunit,'(3e15.7)')st(i),st(im1),st(im2)
                    endif
                
c end of <if (yon_rf.eq.'y')> clause
                else
                    src_stat%yon_ct = 'n'
                endif

c end of <do while (src_stat%yon_ct .eq. 'y')> loop
            enddo

c
c MP3- set output value for stp, std and xyzf
            src_stat%stp = st(im1)
c           write(*,*)'lsearch_p: 5'
            call MPI_BCAST(src_stat, 1, MPI_SEARCH_STATUS,
     >                     0, MPI_COMM_WORLD, ierr)
c
            bcn(1:nb) = bc(1:nb)+src_stat%stp*ds(1:nb)
c
            call cpt_dat_vals_p(shdeg, nb, nd,
     >                          nlocpts, nlocdatpts,
     >                          d2a, dra, dlf, ddlf,
     >                          bcn, ppos, xyzf)
c
            call cptstd_dp(npts, nlocpts,
     >                     cov, jcov, ddat,
     >                     xyzf, std)
c
c
c end of <if (rank .eq. 0)> clause
        else
c
c  Sp wait for starting order from master
            call MPI_BCAST(src_stat, 1, MPI_SEARCH_STATUS,
     >                     0, MPI_COMM_WORLD, ierr)
c
c  SP while asked to do some work
            do while (src_stat%yon_ct .eq. 'y')
c  SP  Receive bcn value from master
c  SP do the work
                bcn(1:nb) = bc(1:nb)+src_stat%stp*ds(1:nb)
c
                call cpt_dat_vals_p(shdeg, nb, nd,
     >                              nlocpts, nlocdatpts,
     >                              d2a, dra, dlf, ddlf,
     >                              bcn, ppos, xyzf)
c
                call cptstd_dp(npts, nlocpts,
     >                         cov, jcov, ddat,
     >                         xyzf, std)

c  SP wait for info on next iteration
                call MPI_BCAST(src_stat, 1, MPI_SEARCH_STATUS,
     >                         0, MPI_COMM_WORLD, ierr)
            enddo

c  SP receive final stp from master & does the final piece of work 
            bcn(1:nb) = bc(1:nb)+src_stat%stp*ds(1:nb)
c
            call cpt_dat_vals_p(shdeg, nb, nd,
     >                          nlocpts, nlocdatpts,
     >                          d2a, dra, dlf, ddlf,
     >                          bcn, ppos, xyzf)
c
            call cptstd_dp(npts, nlocpts,
     >                     cov, jcov, ddat,
     >                     xyzf, std)
        endif
c
        deallocate(bcn)
c
        stp = src_stat%stp
c
        return
        end