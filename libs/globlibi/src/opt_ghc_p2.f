ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         subroutine opt_ghc_p2
c                Vincent Lesur 05/11/2007
c
c     derived from the code opt_pr_p2 19.08.2011
c
c     NOTE: to be used only with diagonal covriance matrix
c
c         Parallel Optimizer using Quasi-Newton methods for linear 
c         problems
c         1) - Use a Conjugate Gradient method with preconditioning
c
c     input:
c         path          path where should be writen outputs
c         itmax(3)      array for Maximum number of iterations
c         npmax         number max of data point with correlated errors
c         nd            space dimension
c         nb            Number of parameters
c         proc_np       number of data+sampling points for all ranks
c         ppos          data point position in ndD + data value
c         BC            Estimate of Base function coefficients
c         dl(3)         control process parameter
c         sub_base_i    Base subroutine to use (see mkArows.f)
c         sub_damp      damping -- not implemented
c         fun_base_f    Base subroutine to use (see cpt_dat_vals_p[2].f)
c         fun_mf        misfit function (like l2_norm.f)
c         fun_std       std function
c         cov(*)        covariance matrix in SLAP Column format
c         jcov          Integer vector describing cov format
c         stdt          target STD value
c
c       output:
c         stdt          STD value for given BC
c         xyzf(*)       Forward modelling for given BC
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine opt_ghc_p2(path, itmax, npmax, nd, nb,
     >                        proc_np, ppos, bc, dl,
     >                        sub_base_i, sub_damp,
     >                        fun_base_f, fun_mf, fun_std,
     >                        cov, jcov, stdt,
     >                        xyzf, bb, gg)
c
        implicit none
c
        include 'mpif.h'
c
        integer itmax(*),npmax,nd,nb
        integer proc_np(*)
        real*8 ppos(*),bc(*),dl(*),cov
        integer jcov(*)
        real*8 stdt,xyzf(*)
        real*8, optional :: bb(:),gg(:)
        character path*100
c
        real*8 fun_base_f, fun_mf, fun_std
        external sub_base_i, sub_damp
        external fun_base_f, fun_mf, fun_std
c
        integer i,ip,it,itm,iunit,ipth,itm_r
        integer ierr,rank,nlocpts
        real*8 stdo,stp,std,epss,dd,cond,dm,beta
        character yon*5
c
c       yon(1:1)    => restart iteration
c       yon(2:2)    => forward modelling
c       yon(3:3)    => gj,ghj calculation
c       yon(4:4)    => GC iteration
c       yon(5:5)    => bc upgrade
c
        real*8, allocatable :: ds(:),dh(:),ddat(:),zz(:)
        real*8, allocatable :: gj(:),gjo(:)
        real*8, allocatable :: ghj(:),ghjo(:)
c
c All defining parallel enviroment
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
        nlocpts = proc_np(rank+1)
c
        if (rank.eq.0) allocate (gjo(1:nb),ghjo(1:nb))
        allocate (zz(1:nb))
        allocate (ds(1:nb))
        allocate (gj(1:nb))
        allocate (ghj(1:nb))
        allocate (dh(1:nb))
        allocate (ddat(1:nlocpts))
c
c Open file for linear search outputs
        if(rank.eq.0) then
            ipth=0
            do while (path(ipth+1:ipth+1).ne.' ')
                ipth=ipth+1
            enddo
            iunit=17
            open(iunit,file=path(1:ipth)//'opt_ghc_p.log')
        endif
c
c Initial
c I assume that BC is the same for all processes at this point
        it=0
        itm=itmax(1)       !total number of iterations
        itm_r=itmax(3)     !number of iterations before restart
        stdo=1.d99
        std=1.d99
c
        epss=dl(1)
        cond=dl(3)
        if (cond.le.0.0d0) cond=1.d30
c
c Starting conditions
        yon(1:5)='yyyyi'
        if (itm.eq.0) yon(3:5)='nnn'
c
        dh(1:nb)=0.0d0
        gj(1:nb)=0.0d0
        ghj(1:nb)=0.0d0
        if (rank.eq.0) then
            gjo(1:nb)=0.0d0
            ghjo(1:nb)=0.0d0
        endif
c
c All define data set
        do ip=1,nlocpts
            ddat(ip)=ppos(ip*(nd+1))
        enddo
c
c All start iteration
        do while (yon(1:1).eq.'y')
            it=it+1
c
c All: do their part in forward modelling
            if (yon(2:2).eq.'y') then
                stdo=std
                call cpt_dat_vals_p2(nd, nlocpts,
     >                               ppos, nb, bc,
     >                               fun_base_f, xyzf)
                call cptstd_dp(npmax, proc_np,
     >                         jcov, cov, ddat, xyzf,
     >                         fun_std, std)
            endif
c
c All: do their part in finding GJ, DH
            if (yon(3:3).eq.'y') then
                ip=1
                call ssqgh_dp(npmax, nd, nlocpts,
     >                        ppos, nb, fun_mf, sub_base_i, bc,
     >                        jcov, cov, ddat, xyzf,
     >                        gj, dh)
c All: check ZEROgradiant
                ip=0
                dm=maxval(dabs(gj(1:nb)))/cond
                do i=1,nb
                    if (dabs(gj(i)).lt.dm) then
                        gj(i)=0.0d0
                        ip=ip+1
                    endif
                enddo

                if (rank.eq.0) write(iunit,*) 
     >              'Number of rejected grad =',ip

c All: update using inverse of "Hessien matrix diag"
                dm=maxval(dabs(dh(1:nb)))/cond
                do i=1,nb
                    dh(i)=dmax1(dh(i),dm)
                    ghj(i)=gj(i)/dh(i)
                enddo
            endif
c
c All: GC iteration
            if (yon(4:4).eq.'y') then
c MP:  Set descent direction
                if (rank.eq.0) then
                    dd=dot_product(gjo,ghjo)
                    if (dabs(dd).gt.epss) then
                        dd=dot_product(gj,ghj)/dd
                        ds(1:nb)=ghj(1:nb)+dd*ds(1:nb)
                    else
                        dd=0.0d0
                        ds(1:nb)=ghj(1:nb)
                    endif
c MP:  Check on direction
                    dd=dsqrt(dot_product(ghj,ghj))
                    dd=dd*dsqrt(dot_product(ds,ds))
                    dd=dot_product(ds,ghj)/dd
                    write(*,*)
     >                  'descent . gradient directions =',dd
                    write(*,*)
     >                  'Current STD before linear search = ',std
                    write(iunit,*)
     >                  'descent . gradient directions =',dd
                    write(iunit,*)
     >                  'Current STD before linear search = ',std
                endif
c
c ALL: Master broadcast ds
                call MPI_BCAST(ds(1:nb), nb, MPI_DOUBLE_PRECISION,
     >                         0, MPI_COMM_WORLD,ierr)
c
c ALL: Find GC step
                stp=0.d0
c ALL: compute zz=2.A^t.W.A.ds
                ip=1
                call cptAtWAds_p(npmax, nd, nlocpts,
     >                           ppos, nb, fun_mf, sub_base_i, bc, ds,
     >                           jcov, cov, ddat,
     >                           xyzf, zz)
c MP:  compute step
                if (rank.eq.0) then
                    stp=dot_product(gj(1:nb),ghj(1:nb))        
                    stp=stp/dot_product(ds(1:nb),zz(1:nb))
c
c MP:  Save cuurent values of gj,ghj
                    gjo(1:nb)=gj(1:nb)
                    ghjo(1:nb)=ghj(1:nb)
c
c MP : update gj, ghj
                    gj(1:nb)=gj(1:nb)-stp*zz(1:nb)
c All: check ZEROgradiant
                    ip=0
                    dm=maxval(dabs(gj(1:nb)))/cond
                    do i=1,nb
                        if (dabs(gj(i)).lt.dm) then
                            gj(i)=0.0d0
                            ip=ip+1
                        endif
                        ghj(i)=gj(i)/dh(i)
                    enddo
                    if (rank.eq.0) write(iunit,*)
     >                  'Number of rejected intermediate grad =',ip
                endif
            endif
c
c MP: What to do next
            if (rank.eq.0) then
                if (yon(1:4).ne.'yynn') then
                    dd=dsqrt(dot_product(gjo,gjo))
                    dd=dd-dsqrt(dot_product(gj,gj))
                    write(*,*)
                    write(*,'(A,I3,x,2e17.9,x,A)')
     >                  'opt_ghc ',it,stp,dd,yon
                    write(iunit,*)
                    write(iunit,'(A,I3,x,2e17.9,x,A)')
     >                  'opt_ghc ',it,stp,dd,yon
c
c MP:  Update if improvement otherwise restart
                    if (dd.gt.0) then
                        yon(1:5)='ynnyy'
                    else
                        gjo(1:nb)=0.0d0
                        ghjo(1:nb)=0.0d0
                        yon(1:5)='yyyyn'
                    endif
c
c MP:  Restart condition
                    if (mod(it,itm_r).eq.0) then
                        gjo(1:nb)=0.0d0
                        ghjo(1:nb)=0.0d0
                        yon(1:4)='yyyy'
                    endif
c
c MP:  Out conditions
c       by max iteration
                    if (it.ge.itm) yon(1:4)='yynn'
c       by no improvement
                    if ((stdo-std)/stdo.lt.epss) yon(1:4)='nnnn'
c       by target achieved
                    if (std.lt.stdt) yon(1:4)='nnnn'
                else
                    yon(1:5)='nnnnn'
                endif
c
                write(*,*)'opt_ghc_p next-step : ',yon
                write(*,*)
                write(iunit,*)'opt_ghc_p next-step : ',yon
                write(iunit,*)
c
c MP: Update parameters
                if (yon(5:5).eq.'y') then
                    bc(1:nb)=bc(1:nb)+stp*ds(1:nb)
                endif
c
c MP: save temp model file
                open(11,file=path(1:ipth)//'model_temp.out')
                    write(11,'(A,e15.7)')'# STD =',std
                    write(11,'(A,i4)')'# IT =',it
                    write(11,'(2A)')'# yon =',yon(1:5)
                    write(11,'(A)')'###########'
                    do i=1,nb
                        write(11,'(i8,2e15.7)')i,bc(i),stp*ds(i)
                    enddo
                close(11)
            endif
c
c ALL: MP Broadcast the information to SPs if required
            call MPI_BCAST(yon, 5, MPI_CHARACTER,
     >                     0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(bc, nb, MPI_DOUBLE_PRECISION,
     >                     0, MPI_COMM_WORLD, ierr)

            if (yon(3:4).eq.'ny') then
                call MPI_BCAST(gj, nb, MPI_DOUBLE_PRECISION,
     >                         0, MPI_COMM_WORLD, ierr)
                call MPI_BCAST(ghj, nb, MPI_DOUBLE_PRECISION,
     >                         0, MPI_COMM_WORLD, ierr)
            endif
        enddo
c
        if (rank.eq.0) close(iunit)
        deallocate (ddat)
        deallocate (dh)
        deallocate (ghj)
        deallocate (gj)
        deallocate (ds)
        deallocate (zz)
        if (rank.eq.0) deallocate (gjo,ghjo)
c
        stdt=std
c        
        return
        end