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
c          2- GG & BB added as defined optional for compatibility with
c             opt_drw_p*
c          3- call to lsearch_p modified
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
c    if itmax(1) < 0 the first descent direction is given by BB
c    if itmax(3) < 0 never choose GC step, uses linear search
c
c     input:
c         path          path where should be writen outputs
c         itmax(3)      array for Maximum number of iterations
c         npmax         number max of data point with correlated errors
c         npm           Maximum total number of data points
c         nd            space dimension
c         npt           Total Number of data
c         ppos          data point position in ndD + data value
c         nn(*)         Number or base function to use
c         BC            Estimate of Base function coefficients
c         dl(3)         control process + damping factor
c         FM            misfit function (like l2_norm.f)
c         BS            Base subroutine to use
c         FSTD          std Function
c         FDAMP         damping -- not implemented
c         nt(*)         data type 
c         cov(*)        covariance matrix in SLAP Column format
c         icov/jcov     Integer vector describing cov format
c         stdt          target STD value
c
c       output:
c         stdt          STD value for given BC
c         xyzf(*)       Forward modelling for given BC
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine opt_pr_p3(path,itmax,npmax,npm,nd,npt,ppos,nn,bc
     >     ,dl,FM,BS,FSTD,FDAMP,nt,cov,icov,jcov,stdt,xyzf,BB,GG)
c
        implicit none
c
        include 'mpif.h'
c
        integer itmax(*),npmax,npm,nd,npt,nn(*),nt(*),icov(*),jcov(*)
        real*8 ppos(*),bc(*),cov(*),stdt,xyzf(*),dl(*)
        character path*100
c
        integer i,ip,it,itm,iunit,ipth,itm_l,itm_r,nb
        real*8 stdo,stp,std,epss,dd,cond,dm,beta
        character yon*5
c
        real*8, optional :: GG(:),BB(:)
c
c       yon(1:1)    => yon_it
c       yon(2:2)    => yon_fwd
c       yon(3:3)    => yon_dwn
c       yon(4:4)    => yon_upg
c       yon(5:5)    => yon_ds: g=GC, p=PR, r=restart, n=restart with stop option
c
        real*8, allocatable :: ds(:),dh(:),ddat(:)
        real*8, allocatable :: gj(:),gjo(:)
        real*8, allocatable :: ghj(:),ghjo(:)
c
        real*8 FM,FSTD
        external FM,BS,FSTD,FDAMP
c
c All defining parallel enviroment
        integer ierr,rank
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
        nb=nn(1)
c
        if(rank.eq.0) then
          allocate (gjo(1:nb))
          allocate (ghjo(1:nb))
        endif
        allocate (ds(1:nb))
        allocate (gj(1:nb))
        allocate (ghj(1:nb))
        allocate (ddat(1:npt))
c
c Open file for linear search outputs
        if(rank.eq.0) then
          ipth=0
          do while (path(ipth+1:ipth+1).ne.' ')
            ipth=ipth+1
          enddo
          iunit=17
          open(iunit,file=path(1:ipth)//'opt_pr_p.log')
        endif
c
c Initial
c I assume that BC is the same for all processes at this point
        it=0
        itm=abs(itmax(1))  !number of iteration
        itm_l=itmax(2)     !number of iteration in linear search
        itm_r=abs(itmax(3))!number of GC iteration between PR
        if(itm_r.le.0)itm_r=1
c
        stdo=1.d99
c
        epss=dl(1)
        cond=dl(3)
        if (cond.le.0.0d0)cond=1.d30
c
c Starting conditions
        yon(1:5)='yyyyr'
        gj(1:nb)=0.0d0
        ghj(1:nb)=0.0d0
        if(rank.eq.0) then
          gjo(1:nb)=0.0d0
          ghjo(1:nb)=0.0d0
        endif
        if(itm.eq.0) yon(3:4)='nn'
c
c All define data set
        do ip=1,npt
          ddat(ip)=ppos(ip*(nd+1))
        enddo
c
c All start iteration
        do while (yon(1:1).eq.'y')
          it=it+1
c All: do their part in forward modelling
          if(yon(2:2).eq.'y')then
            call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,BC,BS,xyzf)
            call cptstd_dp(npmax,npt,nt,icov,jcov,cov,ddat,xyzf
     >                    ,FSTD,std)
            stdo=std
          endif
c
c All: do their part in finding next step length and direction
          if(yon(3:3).eq.'y')then
c All: do their part in finding GJ, DH
            ip=1
            allocate(dh(1:nb))
            dh(1:nb)=0.0d0
c
            if(itmax(1).ge.0.or.it.ne.1)then
              call ssqgh_dp(npmax,npm,npt,ip,nd,ppos,nb,FM,BS,bc
     >               ,icov,jcov,cov,ddat,nt,xyzf,gj,dh)
            else
              if(present(BB))then
                gj(1:nb)=BB(1:nb)
                dh(1:nb)=1.d0
              else
                if(rank.eq.0)
     >          write(*,*)'ERROR: BB required for itmax(1)<0'
                stop
              endif
            endif
c All: check ZEROgradiant
            ip=0
            dm=maxval(dabs(gj(1:nb)))/cond
            do i=1,nb
              if(dabs(gj(i)).lt.dm) then
                 gj(i)=0.0d0
                 ip=ip+1
              endif
            enddo
            if(rank.eq.0)write(iunit,*)'Number of rejected grad =',ip
c All: update using inverse of "Hessien matrix diag"
            dm=maxval(dabs(dh(1:nb)))/cond
            do i=1,nb
              dh(i)=dmax1(dh(i),dm)
              ghj(i)=gj(i)/dh(i)
            enddo
c
            deallocate(dh)
c MP:  Set descent direction
            if(rank.eq.0) then
              beta=0.0d0
              if(yon(5:5).eq.'r'.or.yon(5:5).eq.'n')then
c          SD direction
                ds(1:nb)=ghj(1:nb)
              elseif(yon(5:5).eq.'g')then
c          GC direction
                dd=dot_product(gjo,ghjo)
                if(dabs(dd).gt.epss)then
                  beta=dot_product(gj,ghj)/dd
                  ds(1:nb)=ghj(1:nb)+beta*ds(1:nb)
                else
                  ds(1:nb)=ghj(1:nb)
                endif
              else
c          Polak-Ribiere direction
                dd=dot_product(gjo,ghjo)
                if(dabs(dd).gt.epss)then
                  ghjo(1:nb)=ghj(1:nb)-ghjo(1:nb)
                  beta=dot_product(gj,ghjo)/dd
                  beta=dmax1(beta,0.0d0)
                  ds(1:nb)=ghj(1:nb)+beta*ds(1:nb)
                else
                  ds(1:nb)=ghj(1:nb)
                endif
              endif
c MP:  Check on direction
              dd=dsqrt(dot_product(ghj,ghj))
              dd=dd*dsqrt(dot_product(ds,ds))
              dd=dot_product(ds,ghj)/dd
              write(*,*)'descent . gradient directions =',dd,beta
              write(*,*)'Current STD before linear search = ',std
              write(iunit,*)'descent . gradient directions =',dd,beta
              write(iunit,*)'Current STD before linear search = ',std
            endif
c
c ALL: Master broadcast ds
            call MPI_BCAST(ds(1:nb),nb,MPI_DOUBLE_PRECISION,0
     >                                       ,MPI_COMM_WORLD,ierr)
c ALL: search minimum in descent direction
            stp=0.d0
            if(yon(5:5).eq.'g'.or.yon(5:5).eq.'r')then
              if(itmax(3).ge.0)then
              call gc_step_p(iunit,npm,npmax,nd,npt,ppos,ddat
     >               ,nb,bc,FSTD,BS,nt,cov,icov,jcov,std,gj,ghj,ds
     >               ,stp,xyzf)
              else
              call lsearch_p(iunit,itm_l,npmax,npm,nd,npt,ppos,ddat
     >         ,nb,bc,dl,BS,FSTD,nt,cov,icov,jcov,std,ds
     >         ,stp,xyzf)
              endif
            else
              call lsearch_p(iunit,itm_l,npmax,npm,nd,npt,ppos,ddat
     >         ,nb,bc,dl,BS,FSTD,nt,cov,icov,jcov,std,ds
     >         ,stp,xyzf)
            endif
          endif
c
c MP: What to do next
          if(rank.eq.0)then
          if(yon(3:3).eq.'y')then
            dd=dabs(stp)*dsqrt(dot_product(ds,ds))
            dd=dd/dsqrt(dot_product(bc(1:nb),bc(1:nb)))
            write(*,*)
            write(*,'(A,I3,3e17.9,x,A)')
     >           'opt_pr_p',it,std,dd,(stdo-std)/stdo,yon
            write(iunit,*)
            write(iunit,'(A,I3,3e17.9,x,A)')
     >           'opt_pr_p',it,std,dd,(stdo-std)/stdo,yon
c
c      Out conditions
            if((stdo-std)/stdo.lt.epss)then
              if(yon(5:5).eq.'n')then
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
c
c               if(dabs(stp).gt.epss)then
                  yon(1:5)='yynyn'
                  if(std.gt.stdo)yon(1:5)='yynnn'
c
c      Other conditions
c        Re-start: Failing to find acceptable step
c               else
c                 gjo(1:nb)=0.0d0
c                 ghjo(1:nb)=0.0d0
c                 yon(1:5)='ynyyn'
c                 if(std.ge.stdo)yon(1:5)='yyynn'
c               endif
c        Re-start with stop option: Too small improvement
              else
                yon(1:5)='ynyyn'
                if(std.gt.stdo) then
                  gjo(1:nb)=0.0d0
                  ghjo(1:nb)=0.0d0
                  yon(1:5)='yyynn'
                endif
              endif
            else
c
c        Normal PR or GC iteration
              yon(1:5)='ynyyp'
              if(itmax(3).ge.0) yon(1:5)='ynyyg'
c
c        Re-start: Too small step ... SD iteration
              if(dd.lt.epss)yon(5:5)='r'
c
c        Re-start: No improvement ... restart with stop option
              if(std.gt.stdo)yon(1:5)='yyynn'
c
c        Re-start: SD Iteration
              if(mod(it,itm_r).eq.0)yon(5:5)='r'
            endif
c
c      Emergency stop 
            if(it.ge.itm)yon(1:5)='yynnn'
            if(std.lt.stdt)yon(1:5)='yynnn'
          else
            dd=0.0d0
            yon(1:1)='n'
          endif
          write(*,*)'opt_pr_p next-step : ',yon
          write(*,*)
          write(iunit,*)'opt_pr_p next-step : ',yon
          write(iunit,*)
c
c MP: Update parameters
          if(yon(4:4).eq.'y')then
            bc(1:nb)=bc(1:nb)+stp*ds(1:nb)
            gjo(1:nb)=gj(1:nb)
            ghjo(1:nb)=ghj(1:nb)
            stdo=std
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
c ALL: MP Broadcast the information to SPs
          call MPI_BCAST(yon,5,MPI_CHARACTER,0
     >                                      ,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(bc,nb,MPI_DOUBLE_PRECISION,0
     >                                      ,MPI_COMM_WORLD,ierr)
        enddo
c
        stdt=std
        if(rank.eq.0)close(iunit)
c
        deallocate (ddat)
        deallocate (ghj)
        deallocate (gj)
        deallocate (ds)
        if(rank.eq.0) then
          deallocate (ghjo)
          deallocate (gjo)
        endif
c        
        return
        end
