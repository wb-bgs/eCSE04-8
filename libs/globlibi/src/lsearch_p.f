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
c         npmax         number max of data point with correlated errors
c         npm           Maximum total number of data points
c         nd            space dimension
c         npt           Total Number of data
c         ppos          data point position in ndD
c         ddat          data values
c         nb            Number or base function to use
c         BC            Estimate of Base function coefficients
c         dl(3)         control lsearch process & damping
c         BS            Base subroutine to use
c         FSTD          std Function
c         nt(*)         data type 
c         cov(*)        covariance matrix in SLAP Column format
c         icov/jcov     Integer vector describing cov format
c         std           STD value for given BC
c         ds(*)         gradient direction
c
c       output:
c         stp           recommended step in direction ds(*)
c         std           STD value for given BC+stp*DS
c         xyzf(*)       Forward modelling for given BC+stp*DS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine lsearch_p(iunit,itm,npmax,npm,nd,npt,ppos,ddat
     >  ,nb,bc,dl,BS,FSTD,nt,cov,icov,jcov,std,ds,stp,xyzf)
c
        implicit none
c
        include 'mpif.h'
c
        integer iunit,itm,npmax,npm,nd,npt,nb,nt(*),icov(*),jcov(*)
        real*8 ppos(*),ddat(*),bc(*),cov(*),std,xyzf(*),ds(*),dl(*),stp
c
        integer i,im1,im2,it
        real*8 dj(3),st(3),fct,fctt,fctl,fcts,dd,rt
        real*8 epss,stp1,stp2,gr,gl
        real*8, allocatable :: bcn(:)
        character yon_ct,yon_rf
c
        real*8 FSTD
        external BS,FSTD
c
c  All defining parallel enviroment
        integer ierr,rank
        call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
c
        allocate(bcn(1:nb))
c
c MP: does the all lot
        if(rank.eq.0) then
c
c MP0-Initial
        yon_ct='y'
        epss=dl(1)
        fctl=5.d0
        fcts=1.d0/fctl
        rt=4.d0
c
        fct=2.d0
        if(stp.eq.0.0d0)stp=1.d0
        stp=stp/fct
c
c MP1-Find the large steps
c
c MP 11- initial (given in input: “std” no need to compute)
        dj(1)=std
        st(1)=0.0d0
        dj(2)=std
        st(2)=0.0d0
c
c MP 12- find st(3) such that dj(3) significantly different from dj(1)
        dd=0.d0
        it=1
        do while (dd.le.epss.and.it.le.itm)
          stp=stp*fct
          call MPI_BCAST(yon_ct,1,MPI_CHARACTER,0
     >                                       ,MPI_COMM_WORLD, ierr)
          call MPI_BCAST(stp,1,MPI_DOUBLE_PRECISION,0
     >                                       ,MPI_COMM_WORLD, ierr)
          bcn(1:nb)=bc(1:nb)+stp*ds(1:nb)
          call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,bcn,BS,xyzf)
          call cptstd_dp(npmax,npt,nt,icov,jcov,cov,ddat,xyzf,FSTD,std)
          dj(3)=std
          dd=dabs(dj(3)-dj(1))/dj(1)
          it=it+1
        enddo
        if(it.gt.itm) then
          write(iunit,*)'Lsearch: cannot find minimum step'
          yon_ct='n'
          im1=1
        endif
        st(3)=stp
c
c MP 13- Find st(i),st(im1),st(im2)such that dj(im1)<dj(i) & dj(im1)<dj(im2) 
        if(yon_ct.eq.'y')then
          it=1
          if (dj(3).lt.dj(1)) then
            im1=2
            i=3
            fctt=fct
            do while (dj(i).lt.dj(im1).and.it.le.itm)
              im1=i
              i=mod(im1,3)+1
              stp=stp*fctt
              fctt=fctt*fct
              call MPI_BCAST(yon_ct,1,MPI_CHARACTER,0
     >                                       ,MPI_COMM_WORLD, ierr)
              call MPI_BCAST(stp,1,MPI_DOUBLE_PRECISION,0
     >                                       ,MPI_COMM_WORLD, ierr)
              bcn(1:nb)=bc(1:nb)+stp*ds(1:nb)
              call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,bcn,BS,xyzf)
              call cptstd_dp(npmax,npt,nt,icov,jcov,cov
     >                                          ,ddat,xyzf,FSTD,std)
              dj(i)=std
              st(i)=stp
              it=it+1
            enddo
            im2=6-(i+im1)
            if(it.gt.itm) then
              write(iunit,*)'Lsearch: cannot increase functional value'
              yon_ct='n'
            endif
          else
            im1=3
            im2=1
            fctt=fct
            stp=-stp/fctt
            do while (dj(im2).lt.dj(im1).and.it.le.itm)
              im1=im2
              im2=mod(im1,3)+1
              stp=stp*fctt
              fctt=fctt*fct
              call MPI_BCAST(yon_ct,1,MPI_CHARACTER,0
     >                                       ,MPI_COMM_WORLD, ierr)
              call MPI_BCAST(stp,1,MPI_DOUBLE_PRECISION,0
     >                                       ,MPI_COMM_WORLD, ierr)
              bcn(1:nb)=bc(1:nb)+stp*ds(1:nb)
              call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,bcn,BS,xyzf)
              call cptstd_dp(npmax,npt,nt,icov,jcov,cov
     >                                          ,ddat,xyzf,FSTD,std)
              dj(im2)=std
              st(im2)=stp
              it=it+1
            enddo
            i=6-(im2+im1)
            if(it.gt.itm) then
              write(iunit,*)'Lsearch: cannot increase functional value'
              yon_ct='n'
            endif
          endif
          write(iunit,*)'lsearch start',it
          write(iunit,'(A,3e15.7)')'lsearch start',dj(i),dj(im1),dj(im2)
          write(iunit,'(A,3e15.7)')'lsearch start',st(i),st(im1),st(im2)
        endif
c
c MP2-Find the zero gradient
        it=1
        do while (yon_ct.eq.'y')
          it=it+1
          stp1=st(i)-st(im1)  
          stp2=st(im1)-st(im2)
c
c MP 20- check for step size
          yon_rf='y'
          if(dabs(st(im1)).gt.epss)then
            dd=dabs((stp1+stp2)/st(im1))
          else
            dd=dabs(stp1+stp2)
          endif
          if(dd.lt.epss)yon_rf='n'
c
c MP 21- If step size large enough: Find the next step try
          if(yon_rf.eq.'y')then
c MP  211- check for step interval relative sizes and set default values if needed 
            fct=dabs(stp1/stp2)
            if(fct.lt.fcts) then
              stp=st(im1)-(st(im1)-st(im2))/rt
            elseif(fct.gt.fctl)then
              stp=st(im1)+(st(i)-st(im1))/rt
            else
c MP  212- steps are similar=> computes gradients (gr*gl<0)
c        stp given by False position algorithm
              gr=(dj(i)-dj(im1))/stp1
              gl=(dj(im1)-dj(im2))/stp2
c
              stp=(st(im2)*gr-st(i)*gl)/(gr-gl)
              stp=0.5d0*(st(im1)+stp)
c
c MP  213-check step size and set default values if needed
              if(dabs((stp-st(im1))/max(abs(st(im1)),abs(stp),epss))
     >         .lt.epss) then
                if(fct.lt.1.d0) then
                  stp=st(im1)-(st(im1)-st(im2))/rt
                else
                  stp=st(im1)+(st(i)-st(im1))/rt
                endif
              endif
            endif
c
c MP 22-  Find functional values and new im1,im2 & i
            call MPI_BCAST(yon_ct,1,MPI_CHARACTER,0
     >                                       ,MPI_COMM_WORLD, ierr)
            call MPI_BCAST(stp,1,MPI_DOUBLE_PRECISION,0
     >                                       ,MPI_COMM_WORLD, ierr)
            bcn(1:nb)=bc(1:nb)+stp*ds(1:nb)
            call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,bcn,BS,xyzf)
            call cptstd_dp(npmax,npt,nt,icov,jcov,cov,ddat,xyzf
     >                                                   ,FSTD,std)
            if (stp.gt.st(im1))then
              if(std.ge.dj(im1))then
                dj(i)=std
                st(i)=stp
              else
                dj(im2)=std
                st(im2)=stp
                im1=im2
                im2=6-(i+im1)
              endif
            else
              if(std.ge.dj(im1))then
                dj(im2)=std
                st(im2)=stp
              else
                dj(i)=std
                st(i)=stp
                im1=i
                i=6-(im2+im1)
              endif
            endif
c
c MP 23- Check if maximum iteration 
            if(it.ge.itm) then
              yon_ct='n'
              write(iunit,*)'Lsearch: Maximum iteration',it
              write(iunit,'(3e15.7)')dj(i),dj(im1),dj(im2)
              write(iunit,'(3e15.7)')st(i),st(im1),st(im2)
            endif
c
c MP 24- Check if significant improvement still possible
            dd=max(dj(im2),dj(i))
            dd=dabs((dd-dj(im1))/dj(im1))
            if(dd.lt.epss) then 
              yon_ct='n'
              write(iunit,*)'Lsearch: No hope of improvement',it
              write(iunit,'(3e15.7)')dj(i),dj(im1),dj(im2)
              write(iunit,'(3e15.7)')st(i),st(im1),st(im2)
            endif
          else
            yon_ct='n'
          endif
        enddo
c
c MP3- set output value for stp, std and xyzf
        stp=st(im1)
        call MPI_BCAST(yon_ct,1,MPI_CHARACTER,0
     >                                       ,MPI_COMM_WORLD, ierr)
        call MPI_BCAST(stp,1,MPI_DOUBLE_PRECISION,0
     >                                       ,MPI_COMM_WORLD, ierr)
        bcn(1:nb)=bc(1:nb)+stp*ds(1:nb)
        call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,bcn,BS,xyzf)
        call cptstd_dp(npmax,npt,nt,icov,jcov,cov,ddat,xyzf,FSTD,std)
c
        else
c
c  Sp wait for starting order from master
          call MPI_BCAST(yon_ct,1,MPI_CHARACTER,0
     >                                       ,MPI_COMM_WORLD, ierr)
c
c  SP while asked to do some work
          do while (yon_ct.eq.'y')
c  SP  Receive bcn value from master
            call MPI_BCAST(stp,1,MPI_DOUBLE_PRECISION,0
     >                                       ,MPI_COMM_WORLD, ierr)
c  SP do the work
            bcn(1:nb)=bc(1:nb)+stp*ds(1:nb)
            call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,bcn,BS,xyzf)
            call cptstd_dp(npmax,npt,nt,icov,jcov,cov,ddat,xyzf
     >                                                     ,FSTD,std)
c  SP wait for info on next iteration
            call MPI_BCAST(yon_ct,1,MPI_CHARACTER,0
     >                                       ,MPI_COMM_WORLD, ierr)
          enddo
c  SP receive final stp from master & does the final piece of work 
          call MPI_BCAST(stp,1,MPI_DOUBLE_PRECISION,0
     >                                       ,MPI_COMM_WORLD, ierr)
          bcn(1:nb)=bc(1:nb)+stp*ds(1:nb)
          call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb,bcn,BS,xyzf)
          call cptstd_dp(npmax,npt,nt,icov,jcov,cov,ddat,xyzf,FSTD,std)
        endif
c
        deallocate(bcn)
c
        return
        end
