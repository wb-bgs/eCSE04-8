cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Program MOD_WMAM
c    
c    Build a spherical harmonic model from a regular grid of 
c    total intensity data in Geocentric
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	program mod_wmam
c
        use sph_wmam
c
        implicit none
c
        include 'mpif.h'
c
        integer nbm,npm,nd,npmax,ncm
        parameter (nbm=4000000,npm=6000000,nd=7,npmax=3,ncm=3*npm)
c
        iNtEgEr DEBUG_BHAM, I_DEBUG
        character cc,nom*100
        integer i,j,k,nb(2),npt,nc,nt(npm),ijcov(ncm,2),itmax(3)
        integer ix,iy,nx,ny,il,im
        real*8  ppos(nd+1,npm),bc(nbm),cov(ncm),dw(npm),stdt,dl(3)
        real*8  std,dd, deg_res
c
        integer nub
        real*8 wgh
c
        real*8, allocatable :: GG(:,:),BB(:)
c
        real*8 dsind
        external dsind
c
        real*8 l2_std,l2_norm
        external l2_std,l2_norm
c
        external sph_bi,damp_rien
c
c  Initialize MPI, determine rank
         integer ierr,rank
         call MPI_Init(ierr)
         call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
c  Settings
        call init_sph_wmam(nb(1),bc)
        nb(2)=nb(1)
        if(rank.eq.0)write(*,*)'Number of parameters :',nb(1)
c
c  Read data
        deg_res = 1
        ny=int(1.0/deg_res)*360
        nx=int(1.0/deg_res)*180-1
c
        nom='./Data/wdmam_geocentric.dat'
        if(rank.eq.0)write(*,'(A,A)')'Filename: ',nom
c
        if(rank.eq.0)write(*,*)'Reading file'
        open(10,file=nom,status='old')
          i=0
          do ix=1,nx
            do iy=1,ny
              i=i+1
              read(10,*)ppos(2,i),ppos(1,i),ppos(3,i),ppos(nd+1,i)
              ppos(1,i)=90.0d0-ppos(1,i)
              ppos(4,i)=ryg
            enddo
          enddo
        close(10)
        npt=i
        if(rank.eq.0)write(*,*)'Number of data 1:',npt
c
c  Define the reference field
        if(rank.eq.0)write(*,*)"define reference field values"
        nom='./Data/coef_1990_15.dat'
        i=nbm
        call read_model(nom,ryg,bc,i)
c
        nb(1)=15*17
c
        nt(1:npt)=1
        dw=0.0d0
        call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb(1),bc,sph_bi,dw)
        ppos(5,1:npt)=dw(1:npt)
        if(rank.eq.0)write(*,'(A)')"X CM4 component calculated"
c
        nt(1:npt)=2
        dw=0.0d0
        call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb(1),bc,sph_bi,dw)
        ppos(6,1:npt)=dw(1:npt)
        if(rank.eq.0)write(*,'(A)')"Y CM4 component calculated"
        nt(1:npt)=3
        dw=0.0d0
        call cpt_dat_vals_p(nd,npm,npt,nt,ppos,nb(1),bc,sph_bi,dw)
        ppos(7,1:npt)=dw(1:npt)
        if(rank.eq.0)write(*,'(A)')"Z CM4 component calculated"
c
c  Define covariance matrix: sin(colat) weight
        if(rank.eq.0)write(*,'(A)')"define cov matrix"
        nt(1:npt)=1
        nc=npt
c
        do i=1,nc
          cov(i)=dsind(ppos(1,i))
          ijcov(i,1)=i
          ijcov(i,2)=i
          cov(i)=1.d0/cov(i)
        enddo
c
c  Add smoothing equations
        if(rank.eq.0)write(*,'(A)')"define regularisation"
        nub=100
        wgh=5.0d0
        call build_damp_space(nub,npm,npt,ncm,nc,nd
     >                     ,ilg,wgh,nt,ijcov,cov,ppos)
c
c  finalise covariance matrix
        call DS2Y(npt,npt,ijcov(1,1),ijcov(1,2),cov,0)
c
c  Define starting model
        nb(1)=ilg*(ilg+2)
        nb(2)=nb(1)
c
        if(rank.eq.0)write(*,'(A)')"define starting model"
        bc(1:nb(1))=0.0d0
        open(10,file='./Data/model.in',status='old')
          read(10,*)nom
          read(10,*)nom
          read(10,*)nom
          read(10,*)nom
          do i=1,nb(1)
            read(10,*)j,bc(i),dd
          enddo
        close(10)
c
c  Invert data
        itmax(1)=7
        itmax(2)=5
        itmax(3)=-10
        stdt=1.0d0
        dl(1)=1.0d-10
        dl(2)=0.0d0
        dl(3)=1.d14
c
        dw=1.d0
c
c       call flush()
c
        if(rank.eq.0) then
          write(*,'(A)')"Start Inversion"
          write(*,*)"ilg : ",ilg
          write(*,*)"nb : ",nb(1),nb(2)
          write(*,*)"npt : ",npt
          write(*,*)"wgh : ",wgh
          write(*,*)"Itmax : ",itmax(1:3)
          write(*,*)"dl : ",dl(1:3)
        endif
c
        nom='./Results/'
        allocate(GG(1:1,1:1))
        allocate(BB(1:1))
c
        call opt_pr_p3(nom,itmax,npmax,npm,nd,npt,ppos,nb,bc,dl
     >          ,l2_norm,sub_sph_wmam_l,l2_std,damp_rien
     >          ,nt,cov,ijcov(1,1),ijcov(1,2),stdt,dw,BB,GG)
c       call opt_ghc_p2(nom,itmax,npmax,npm,nd,npt,ppos,nb,bc,dl
c     >          ,l2_norm,sub_sph_wmam_l,l2_std,damp_rien
c     >          ,nt,cov,ijcov(1,1),ijcov(1,2),stdt,dw,BB,GG)

c       call opt_pr_p3(nom,itmax,npmax,npm,nd,npt,ppos,nb,bc,dl
c    >          ,wmam_norm,sub_sph_wmam_l,wmam_var,damp_rien
c    >          ,nt,cov,ijcov(1,1),ijcov(1,2),stdt,dw,BB,GG)
c
        deallocate(BB)
        deallocate(GG)
c
        if(rank.eq.0)then
          write(*,'(A)')' '
          write(*,'(A,e15.7)')'The L2 STD is: ',stdt
          write(*,'(A)')' '
c
c  Writing fit to data per component
          nom='./Results/fit_No_P.out'
          call write_comp(nom,npt,npm,nd,nt,1,ppos,dw)
          nom='./Results/fit_damp.out'
          call write_comp(nom,npt,npm,nd,nt,100,ppos,dw)
c
          std=stdt*npt/(npt-nb(1))
          dw(1:nb(1))=std
c
c  Saving update base coefficients
          nom='./Results/model_No_P.out'
          open(10,file=nom)
            write(10,'(A)')'#'
            write(10,'(A,i8)')'#lmax= ',ilg
            write(10,'(A)')'###########'
            write(10,*)ryg
            i=0
            do il=1,ilg
              im=0
              i=i+1
              write(10,*)'c',il,im,bc(i),dw(i)
              do im=1,il
                i=i+1
                write(10,*)'c',il,im,bc(i),dw(i)
                i=i+1
                write(10,*)'c',il,-im,bc(i),dw(i)
              enddo
            enddo
          close(10)
        endif
c
        call MPI_Finalize(ierr)
c
        stop
        end
