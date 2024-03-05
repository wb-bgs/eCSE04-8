ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       subroutine build_damp.f   V.Lesur
c
c       Setup equations for imposing the smoothing Norm:
c        int_{\Omega} \nabla B_r \dot \nabla B_r d_{\omega}
c
c       input:
c       nub       integer    index for this kind of data
c       npm       integer    maximum number of data points acceptable
c       npt       integer    current number of data points available
c       ncm       integer    maximum number of covarience values
c       nc        integer    current number of covarience values
c       nd        integer    such that nd+1 is the lead dim of ddat
c       llm       integer    maxium degree for lithospheric field
c       wgh       real*8     weight for damping
c
c       output:
c       npt       integer    current number of data points available
c       nc        integer    current number of covariance values
c       nt(*)     integer    array seting data type to nub
c       ddat(*)   real*8     array of data values
c       ijcov(*)  integer    array defining cov matrix
c       cov(*)    real*8     array of cov matrix
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine build_damp(nub,npm,npt,ncm,nc,nd,
     >                        llm,wgh,nt,ijcov,cov,ddat)
c
        use sph_wmam
c
        implicit none
c
        real*8, parameter :: RAG=6371.2d0
c
        integer nub,npm,npt,ncm,nc,nd,llm,ijcov(ncm,*),nt(*)
        real*8 ddat(nd+1,*),cov(*),wgh
c
        integer ll,ml
c
c loop over the degrees & orders
        do ll=1,llm
          ml=0
          npt=npt+1
          ddat(1:nd+1,npt)=0.0d0
          ddat(1,npt)=dble(ll)
          ddat(2,npt)=dble(ml)
          ddat(3,npt)=RAG
          nt(npt)=nub
c
          nc=nc+1
          ijcov(nc,1)=npt
          ijcov(nc,2)=npt
          cov(nc)=1.d0/dble((ll+1)*(ll+1))
          cov(nc)=cov(nc)*dble(2*ll+1)/dble(ll*(ll+1))
          cov(nc)=cov(nc)/wgh
c
          do ml=1,ll
            npt=npt+1
            ddat(1:nd+1,npt)=0.0d0
            ddat(1,npt)=dble(ll)
            ddat(2,npt)=dble(ml)
            ddat(3,npt)=RAG
            nt(npt)=nub
c
            nc=nc+1
            ijcov(nc,1)=npt
            ijcov(nc,2)=npt
            cov(nc)=1.d0/dble((ll+1)*(ll+1))
            cov(nc)=cov(nc)*dble(2*ll+1)/dble(ll*(ll+1))
            cov(nc)=cov(nc)/wgh
c
            npt=npt+1
            ddat(1:nd+1,npt)=0.0d0
            ddat(1,npt)=dble(ll)
            ddat(2,npt)=dble(-ml)
            ddat(3,npt)=RAG
            nt(npt)=nub
c
            nc=nc+1
            ijcov(nc,1)=npt
            ijcov(nc,2)=npt
            cov(nc)=1.d0/dble((ll+1)*(ll+1))
            cov(nc)=cov(nc)*dble(2*ll+1)/dble(ll*(ll+1))
            cov(nc)=cov(nc)/wgh
          enddo
        enddo
c
        return
        end subroutine build_damp
