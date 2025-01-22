ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         subroutine opt_pr_p3
c                Vincent Lesur 05/11/2007
c
c     NOTE: to be used only with diagonal covriance matrix
c
c  THe switch from GC to PR is not working as expected
c       itmax 3>0  => GC
c       itmax 3<0  => PR
c
c         modified 10.05.2013 V. Lesur
c          1- Cleanup / reorganize allocation step
c         modified 01.11.2012 V.Lesur
c          1- Change choice for PR descent direction
c          1- ITMAX(3)<0 no GC steps are used
c         modified 02.12.2011 V.Lesur
c          1- to give the possibility to overwrite the first descent direction
c         modified 14.07.2011 V.Lesur
c          1- nb defined as an array for compatibility woth opt_drw_p*
c          2- call to lsearch_p modified
c         modified 01.06.2011 V.Lesur
c          1- no out conditions fo stp too small
c          2- counter for rejected gradient directions
c         modified 10.03.2011 V.Lesur
c          1- sign of gj changed
c         modified 16/02/2011 v.Lesur
c          1- Introduction of itmax as an array
c         modified 01/11/2010 (v.lesur) 
c          1- Revise the iteration tests
c          2- Introduce the GC Beta and alpha
c          3- Replace the lsearch_gh_p.f subroutine call by a call
c             to gc_step_p.f 
c
c         modified 05/11/2007 from opt_gh_p.f
c
c         Parallel Optimizer using Quasi-Newton methods
c         1) - Use a Conjugate Gradient method with preconditioning
c         2) - beta calculated following Polak-Ribiere formula
c
c       itmax(1)        overall number of iteration
c       itmax(2)        number of iteration in linear search
c       itmax(3)        number of iteration before restart
c    if itmax(3) < 0 never choose GC step, uses linear search
c
c     input:
c         path          path where should be writen outputs
c         itmax(3)      array for Maximum number of iterations
c         shdeg         max SH degree value
c         nb            Number of parameters
c         nd            space dimension
c         npts          Total number of points (data + sampling) for all ranks
c         nlocpts       Total number of points for this rank
c         nlocdatpts    number of data points assigned to rank
c         bc            Estimate of Base function coefficients
c         ppos          data point position in ndD + data value
c         dl(3)         control process + damping factor
c         cov(*)        covariance matrix in SLAP Column format
c         jcov          Integer vector describing cov format
c         stdt          target STD value
c
c       output:
c         stdt          STD value for given BC
c         xyzf(*)       Forward modelling for given BC
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine opt_pr_p3(path, itmax, shdeg, nb, nd,
     >                       npts, nlocpts, nlocdatpts,
     >                       bc, ppos, dl, cov, jcov,
     >                       stdt, xyzf)
c
        implicit none
c
        include 'mpif.h'
        include 'mpi_status_types.h'
c
        character path*100 
        integer itmax(1:3), npts
c
        integer shdeg, nb, nd
        integer nlocpts, nlocdatpts
c
        integer jcov(1:nlocpts+2)
        real*8 cov(1:nlocpts)
        real*8 bc(1:nb)
        real*8 ppos(1:nd+1,1:nlocpts)
        real*8 xyzf(1:nlocpts)
c
        real*8 stdt, dl(1:3)
c
        integer i, j, k
        integer ip, it, itm, iunit
        integer ipth, itm_l, itm_r
        integer ierr, rank
        real*8 std0, stp, std
        real*8 epss, dd, cond, dm, beta
c
        type(inversion_status) inv_stat
        type(search_status) src_stat
        integer MPI_INVERSION_STATUS, MPI_SEARCH_STATUS
c
c       inv_stat%yon(1:1) => yon_it
c       inv_stat%yon(2:2) => yon_fwd
c       inv_stat%yon(3:3) => yon_dwn
c       inv_stat%yon(4:4) => yon_upg
c       inv_stat%yon(5:5) => yon_ds: g=GC, p=PR, r=restart, n=restart with stop option
c
c  variables for populating d2a array 
        integer nm
        real*8  dnm, d1, d2
c
c  declare arrays used within "XYZsph_bi0.f"
        real*8, allocatable :: d2a(:)
c
c  declare arrays used to handle output from ssqgh_dp() subroutine
        real*8, allocatable :: dh(:)
        real*8, allocatable :: gj(:)
        real*8, allocatable :: ddat(:)
        real*8, allocatable :: ds(:)
        real*8, allocatable :: gj0(:)
        real*8, allocatable :: ghj(:)
        real*8, allocatable :: ghj0(:)
c
c  slatec subroutine
        real*8 dgamln
        external dgamln 
c
c
c All defining parallel enviroment
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
c
c Commit the inversion and gradient search MPI status types
        call init_mpi_status_types(nb,bc,inv_stat,src_stat,
     >                             MPI_INVERSION_STATUS,
     >                             MPI_SEARCH_STATUS)
c
c
c Open file for linear search outputs
        if (rank .eq. 0) then
            ipth = 0
            do while (path(ipth+1:ipth+1) .ne. ' ')
                ipth = ipth+1
            enddo
            iunit = 17
            open(iunit, file=path(1:ipth)//'opt_pr_p.log')
        endif
c
c Initial
c I assume that BC is the same for all processes at this point
        it = 0
        itm = abs(itmax(1))  !number of iteration
        itm_l = itmax(2)     !number of iteration in linear search
        itm_r = abs(itmax(3))!number of GC iteration between PR
        if (itm_r .le. 0) itm_r = 1
c
        std0 = 1.d99
c
        epss = dl(1)
        cond = dl(3)
        if (cond .le. 0.0d0) cond = 1.d30
c
c Starting conditions
        inv_stat%yon(1:5) = 'yyyyr'
        if (itm .eq. 0) inv_stat%yon(3:4) = 'nn'
c
c  Allocate and initialise arrays used to handle output from ssqgh_dp() subroutine
        allocate(ddat(1:nlocpts))
        do ip = 1,nlocpts
          ddat(ip) = ppos(nd+1,ip)
        enddo
c
        allocate(ds(1:nb))
        allocate(gj(1:nb))
        allocate(ghj(1:nb))
        gj(1:nb)  = 0.0d0
        ghj(1:nb) = 0.0d0
        if (rank .eq. 0) then
            allocate(gj0(1:nb))
            allocate(ghj0(1:nb))
            gj0(1:nb)  = 0.0d0
            ghj0(1:nb) = 0.0d0
        endif
c
c  Allocate private arrays used within "XYZshp_bi0.f"
        allocate(d2a(0:shdeg))
c 
c  Initialize d2a array used within mk_lf_dlf() subroutine
        do nm = 0,shdeg
          dnm = dble(nm)
          d1 = dgamln(2*dnm+1.0d0, ierr)
          d2 = dgamln(dnm+1.0d0, ierr)
          if (ierr .ne. 0) then
            write(*,*) 'd2a: Cannot computes normalisation cst !'
            stop
          endif
c
          d2 = 0.5d0*d1 - d2
          d2 = d2 - nm*dlog(2.0d0)
          d2 = dexp(d2)
          if (nm .ne. 0) d2 = d2*dsqrt(2.0d0)
c
          d2a(nm) = d2
        enddo
c
c
c All start iteration
        do while (inv_stat%yon(1:1) .eq. 'y')
            it = it+1
c All: do their part in forward modelling
            if (inv_stat%yon(2:2) .eq. 'y') then
c               if(rank.eq.0)write(*,*)'opt_pr_p3: 1'           
                call cpt_dat_vals_p(shdeg, nb, nd, nlocpts, nlocdatpts,
     >                              d2a, inv_stat%bc, ppos, xyzf)
c
                call cptstd_dp(npts, nlocpts, cov, jcov,
     >                         ddat, xyzf, std)
                std0 = std
            endif
c
c All: do their part in finding next step length and direction
            if (inv_stat%yon(3:3) .eq. 'y') then
c All: do their part in finding GJ, DH
                ip = 1
                allocate(dh(1:nb))
                dh(1:nb) = 0.0d0
c
                if ((itmax(1) .ge. 0) .or. (it .ne. 1)) then
c                   if(rank.eq.0)write(*,*)'opt_pr_p3: 2'
                    call ssqgh_dp(shdeg, nb, nd, nlocpts, nlocdatpts,
     >                            d2a, inv_stat%bc, ppos,
     >                            cov, jcov, ddat, xyzf, gj, dh)
                else
                    if (rank .eq. 0) then
                        write(*,*)
     >                      'ERROR: bb required for itmax(1)<0'
                    endif
                    stop
                endif

c All: check ZEROgradiant
                ip = 0
                dm = maxval(dabs(gj(1:nb)))/cond
                do i = 1,nb
                    if (dabs(gj(i)) .lt. dm) then
                        gj(i) = 0.0d0
                        ip = ip+1
                    endif
                enddo

                if (rank .eq. 0)
     >              write(iunit,*)'Number of rejected grad =',ip
c All: update using inverse of "Hessien matrix diag"
                dm = maxval(dabs(dh(1:nb)))/cond
                do i = 1,nb
                    dh(i) = dmax1(dh(i),dm)
                    ghj(i) = gj(i)/dh(i)
                enddo
c
                deallocate(dh)
c MP:  Set descent direction
                if (rank .eq. 0) then
                    beta = 0.0d0
                    if (inv_stat%yon(5:5) .eq. 'r' .or. 
     >                  inv_stat%yon(5:5) .eq. 'n') then
c          SD direction
                        ds(1:nb) = ghj(1:nb)
                    elseif (inv_stat%yon(5:5) .eq. 'g') then
c          GC direction
                        dd = dot_product(gj0, ghj0)
                        if (dabs(dd) .gt. epss) then
                            beta = dot_product(gj,ghj)/dd
                            ds(1:nb) = ghj(1:nb)+beta*ds(1:nb)
                        else
                            ds(1:nb) = ghj(1:nb)
                        endif
                    else
c          Polak-Ribiere direction
                        dd = dot_product(gj0,ghj0)
                        if (dabs(dd) .gt. epss) then
                            ghj0(1:nb) = ghj(1:nb)-ghj0(1:nb)
                            beta = dot_product(gj,ghj0)/dd
                            beta = dmax1(beta,0.0d0)
                            ds(1:nb) = ghj(1:nb)+beta*ds(1:nb)
                        else
                            ds(1:nb) = ghj(1:nb)
                        endif
                    endif
c MP:  Check on direction
                    dd = dsqrt(dot_product(ghj,ghj))
                    dd = dd*dsqrt(dot_product(ds,ds))
                    dd = dot_product(ds,ghj)/dd
                    write(*,*)
     >                  'descent . gradient directions =',dd,beta
                    write(*,*)
     >                  'Current STD before linear search = ',std
                    write(iunit,*)
     >                  'descent . gradient directions =',dd,beta
                    write(iunit,*)
     >                  'Current STD before linear search = ',std

c end of <if (rank.eq.0)> clause
                endif
c
c ALL: Master broadcast ds
c               if(rank.eq.0)write(*,*)'opt_pr_p3: 3'
                call MPI_BCAST(ds, nb, MPI_DOUBLE_PRECISION,
     >                         0, MPI_COMM_WORLD, ierr)
c
c ALL: search minimum in descent direction
                stp = 0.d0
                if (inv_stat%yon(5:5) .eq. 'g' .or. 
     >              inv_stat%yon(5:5) .eq. 'r') then
                    if (itmax(3) .ge. 0) then
c                       if(rank.eq.0)write(*,*)'opt_pr_p3: 4'
                        call gc_step_p(iunit, shdeg, nb, nd,
     >                                 npts, nlocpts, nlocdatpts,
     >                                 d2a, inv_stat%bc,
     >                                 ppos, ddat, cov, jcov, std,
     >                                 gj, ghj, ds, stp, xyzf)
                    else
c                       if(rank.eq.0)write(*,*)'opt_pr_p3: 5'
                        call lsearch_p(iunit, itm_l, shdeg, nb, nd, 
     >                                 npts, nlocpts, nlocdatpts,
     >                                 d2a, inv_stat%bc,
     >                                 ppos, ddat, src_stat,
     >                                 MPI_SEARCH_STATUS,
     >                                 dl, cov, jcov,
     >                                 std, ds, stp, xyzf)
                    endif
                else
c                   if(rank.eq.0)write(*,*)'opt_pr_p3: 6'
                    call lsearch_p(iunit, itm_l, shdeg, nb, nd,
     >                             npts, nlocpts, nlocdatpts,
     >                             d2a, inv_stat%bc,
     >                             ppos, ddat, src_stat,
     >                             MPI_SEARCH_STATUS,
     >                             dl, cov, jcov,
     >                             std, ds, stp, xyzf)
                endif

c end of <if (inv_stat%yon(3:3).eq.'y')> clause
            endif
c
c MP: What to do next
            if (rank .eq. 0) then
                if (inv_stat%yon(3:3) .eq. 'y') then
                    dd = dabs(stp)*dsqrt(dot_product(ds,ds))
                    dd = dd/dsqrt(dot_product(inv_stat%bc(1:nb),
     >                                        inv_stat%bc(1:nb)))
                    write(*,*)
                    write(*,'(A,I3,3e17.9,x,A)')
     >           'opt_pr_p',it,std,dd,(std0-std)/std0,inv_stat%yon
                    write(iunit,*)
                    write(iunit,'(A,I3,3e17.9,x,A)')
     >           'opt_pr_p',it,std,dd,(std0-std)/std0,inv_stat%yon
c
c      Out conditions
                    if ((std0-std)/std0 .lt. epss) then
                        if (inv_stat%yon(5:5) .eq. 'n') then
c
                            write(iunit,'(A)')'OUT CONDITIONS'
                            write(iunit,'(A,e17.9)')
     >               'STEP : ',stp
                            write(iunit,'(A,e17.9)')
     >               '||ds|| : ',dsqrt(dot_product(ds,ds))
                            write(iunit,'(A,e17.9)')
     >               '||gj|| : ',dsqrt(dot_product(gj,gj))
                            write(iunit,'(A,e17.9)')
     >               '||ghj|| : ',dsqrt(dot_product(ghj,ghj))
                            inv_stat%yon(1:5) = 'yynyn'
                            if (std .gt. std0) then
                                inv_stat%yon(1:5) = 'yynnn'
                            endif
c        Re-start with stop option: Too small improvement
                        else
                            inv_stat%yon(1:5) = 'ynyyn'
                            if (std .gt. std0) then
                                gj0(1:nb) = 0.0d0
                                ghj0(1:nb) = 0.0d0
                                inv_stat%yon(1:5) = 'yyynn'
                            endif
                        endif

                    else
c
c        Normal PR or GC iteration
                        inv_stat%yon(1:5) = 'ynyyp'
                        if (itmax(3) .ge. 0) then
                            inv_stat%yon(1:5) = 'ynyyg'
                        endif
c
c        Re-start: Too small step ... SD iteration
                        if (dd .lt. epss) then
                            inv_stat%yon(5:5) = 'r'
                        endif
c
c        Re-start: No improvement ... restart with stop option
                        if (std .gt. std0) then
                            inv_stat%yon(1:5)='yyynn'
                        endif
c
c        Re-start: SD Iteration
                        if (mod(it,itm_r) .eq. 0) then
                            inv_stat%yon(5:5) = 'r'
                        endif
                    endif
c
c        Emergency stop 
                    if (it .ge. itm) inv_stat%yon(1:5) = 'yynnn'
                    if (std .lt. stdt) inv_stat%yon(1:5) = 'yynnn'

c end of <if (inv_stat%yon(3:3).eq.'y')> clause
                else
                    dd = 0.0d0
                    inv_stat%yon(1:1) = 'n'
                endif
                write(*,*)'opt_pr_p next-step : ',inv_stat%yon
                write(*,*)
                write(iunit,*)'opt_pr_p next-step : ',inv_stat%yon
                write(iunit,*)
c
c MP: Update parameters
                if (inv_stat%yon(4:4) .eq. 'y') then
                    inv_stat%bc(1:nb) = inv_stat%bc(1:nb)+stp*ds(1:nb)
                    gj0(1:nb) = gj(1:nb)
                    ghj0(1:nb) = ghj(1:nb)
                    std0 = std
                endif
c
c MP: save temp model file
                open(11, file=path(1:ipth)//'model_temp.out')
                    write(11,'(A,e15.7)')'# STD =',std
                    write(11,'(A,i4)')'# IT =',it
                    write(11,'(2A)')'# yon =',inv_stat%yon(1:5)
                    write(11,'(A)')'###########'
                    do i = 1,nb
                        write(11, '(i8,2e15.7)')i,
     >                        inv_stat%bc(i), stp*ds(i)
                    enddo
                close(11)

c end of <if (rank.eq.0)> clause
            endif
c
c ALL: MP Broadcast the information to SPs
c           if(rank.eq.0)write(*,*)'opt_pr_p3: 7'
            call MPI_BCAST(inv_stat, 1, MPI_INVERSION_STATUS,
     >                     0, MPI_COMM_WORLD, ierr)
        
c end of <do while (inv_stat%yon(1:1).eq.'y')> loop
        enddo
c
c
        stdt = std
        if (rank .eq. 0) close(iunit)
c
        deallocate(d2a)
c
        deallocate(ddat, ds, gj, ghj)
        if (rank .eq. 0) then
            deallocate(gj0, ghj0)
        endif
c        
c Free the inversion and gradient search MPI status types
        call free_mpi_status_types(nb,bc,inv_stat,
     >                             MPI_INVERSION_STATUS,
     >                             MPI_SEARCH_STATUS)
c
c
        return
        end