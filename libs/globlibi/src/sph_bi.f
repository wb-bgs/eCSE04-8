cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    subroutine sph_bi
c
c	V. Lesur  12/01/2005
c
c   Recursive version of sph_bi  should be compiled with Fortran 90 (or 95)
c
c   for a given base number:
c      1->  X component
c      2->  Y component
c      3->  Z component
c      4->  F component
c      53->  Z component (GEO)
c      100->  scaled potential
c
c   Computes a vector of the "nb" base elements for a given
c   set of base parameters
c
c   input:
c     iof     inverse or forward (character)
c     nub     base number
c     nb0     number of base element from ref model
c     nb      number of base element based on max spherical deg
c     nd      space dimension
c     bp(*)   base parameters
c     bc(*)   base coefficients
c
c   output:
c     be(*)   nb base elements dim min nb
c             be(1) = be.bc
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine sph_bi(iof, nub, nd, nb0, nb, bc, bp, be)
c
        implicit none
c
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
        character iof
        integer nub, nb0, nb, nd
        real*8 bc(nb), bp(nd), be(nb0)
c
        integer nlie
c
c
        nlie = int(sqrt(float(nb0+1))+0.4)-1
        be(1:nb0) = 0.0d0
c
        if (nub .eq. 1) then
            call Xsph_bi(nd, nb0, nlie, bp, be)
        elseif (nub .eq. 2) then
            call Ysph_bi(nd, nb0, nlie, bp, be)
        elseif (nub .eq. 3) then
            call Zsph_bi(nd, nb0, nlie, bp, be)
        else
            write(*,*)'sph_bi:'
            write(*,*)'Only XYZ base functions available'
            stop
        endif
c
        if (iof .eq. 'f') then
            be(1) = dot_product(be(1:nb0), bc(1:nb0))
        endif
c
        return
        end