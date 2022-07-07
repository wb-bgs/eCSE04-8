cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       subroutine set_FG_sampling
c                                    V. Lesur 07/07/2006
c
c       Calculate the sampling points and weights for the
c       Fourier-Gauss sampling theorem on the Sphere
c
c       The number of sampling point is np=(lmax+1)*(2*lmax+1)
c       input:
c          lmax     INTEGER    Maximum spherical harmonic degree
c       output:
c          np       INTEGER    Number of sampling points
c          vrt(2,*) REAL*8     Array of vertex positions
c          wght(*)  Real*8     Weight associated with each vertex
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine set_FG_sampling(lmax,np,vrt,wght)
c
        implicit none
c
        real*8 vrt(2,*),wght(*)
        integer lmax,np
c
        real*8 dw1,dw2
        integer nn,ilo,ila
        real*8, allocatable :: xx(:),ww(:)
c
        real*8 dacosd
        external dacosd
c
        nn=lmax+1
        dw1=-1.0d0
        dw2=1.0d0
c
c  Find the zeros and weight
        allocate (xx(1:nn),ww(1:nn))
        call gauleg(dw1,dw2,xx,ww,nn)
c
c  Defind sampling in Longitude
        dw2=360.0d0/dble(2*nn-1)
c
c  Define the grid
        np=0
        do ila=1,nn
          do ilo=1,2*nn-1
            np=np+1
c
            vrt(1,np)=dacosd(xx(ila))
            vrt(2,np)=dble(ilo-1)*dw2
            wght(np)=ww(ila)
c
          enddo
        enddo
c
        deallocate (xx,ww)
        return
        end
