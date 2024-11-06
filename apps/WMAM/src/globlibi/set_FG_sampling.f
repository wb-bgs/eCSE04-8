cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       subroutine set_FG_sampling
c                                    V. Lesur 07/07/2006
c
c       Calculate the sampling points and weights for the
c       Fourier-Gauss sampling theorem on the Sphere
c
c       The number of sampling point is np=(lmax+1)*(2*lmax+1)
c       input:
c          lmax           INTEGER  Maximum spherical harmonic degree
c          nlocsampts     INTEGER  Number of sampling points local to rank
c          imin_locsampts INTEGER  First index for sample points
c          imax_locsampts INTEGER  Last index for sample points
c       output:
c          vrt(2,*) REAL*8     Array of vertex positions
c          wght(*)  Real*8     Weight associated with each vertex
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine set_FG_sampling(lmax, nlocsampts,
     >                             imin_locsampts, imax_locsampts,
     >                             vrt, wght)
c
        implicit none
c
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
        integer lmax, nlocsampts
        integer imin_locsampts, imax_locsampts
        real*8 vrt(1:2,1:nlocsampts), wght(1:nlocsampts)
c
        real*8 dw1, dw2
        integer nn, ilo, ila, j, np
        real*8, allocatable :: xx(:), ww(:)
c
c
        nn = lmax+1
        dw1 = -1.0d0
        dw2 = 1.0d0
c
c  Find the zeros and weight
        allocate(xx(1:nn), ww(1:nn))
c
        call gauleg(dw1, dw2, xx, ww, nn)
c
c  Defind sampling in Longitude
        dw2 = 360.0d0/dble(2*nn-1)
c
c  Define the grid
        np = 1
        j = 1
        do ila = 1,nn
          do ilo = 1,2*nn-1
            if (np .ge. imin_locsampts .and.
     >          np .le. imax_locsampts) then
              vrt(1,j) = dacos(xx(ila))/D2R
              vrt(2,j) = dble(ilo-1)*dw2
              wght(j) = ww(ila)
              j = j+1
            endif
            np = np+1
          enddo
        enddo
c
        deallocate(xx, ww)
        return
        end