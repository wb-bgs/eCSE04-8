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
c     nb      nomber of base element
c     bp(*)   base parameters
c     bc(*)   base coefficients
c
c   output:
c     be(*)   nb base elements dim min nb
c             be(1) = be.bc
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        recursive subroutine sph_bi(iof,nub,nb,bc,bp,be)
c
        implicit none
c
        real*8, parameter :: D2R = 4.d0*datan(1.d0)/180.d0
c
        character iof
        integer nub,nb,nli,nti,i
        real*8 bp(*),be(*),bc(*),ry,df(3),dnm
        real*8, allocatable :: dw1(:)
c
        real*8 p1
c
c
        p1 = bp(1)*D2R
c
        nli=int(sqrt(float(nb+1))+0.4)-1
        nti=1
        ry=0.0d0
        do i=1,nb
          be(i)=0.0d0
        enddo
c
        if(nub.eq.1)then
          call Xsph_bi(ry,1,nli,nti,bp,be)
        elseif(nub.eq.2)then
          call Ysph_bi(ry,1,nli,nti,bp,be)
        elseif(nub.eq.3)then
          call Zsph_bi(ry,1,nli,nti,bp,be)
        elseif(nub.eq.4)then
          allocate(dw1(1:nb))
          df(1:3)=0.0d0
          do i=1,3
            call sph_bi('i',i,nb,bc,bp,dw1)
            df(i)=dot_product(dw1,bc(1:nb))
            be(1:nb)=be(1:nb)+dw1*df(i)
          enddo
          deallocate(dw1)
          dnm=dot_product(df,df)
          dnm=dsqrt(dnm)
          be(1:nb)=be(1:nb)/dnm
        elseif(nub.eq.53)then
          allocate(dw1(1:nb))
            call Xsph_bi(ry,1,nli,nti,bp,dw1)
            call Zsph_bi(ry,1,nli,nti,bp,be)
            do i=1,nb
              be(i)=be(i)*dcos(p1)
              be(i)=dw1(i)*dsin(p1)-be(i)
            enddo
          deallocate(dw1)
        elseif(nub.eq.100)then
          call Psph_bi(ry,1,nli,nti,bp,be)
        else
          write(*,*)'sph_bi:'
          write(*,*)'Only XYZF & P base functions available'
          stop
        endif
c
        if (iof.eq.'f') then
          be(1) = dot_product(be(1:nb), bc(1:nb))
        endif
c
        return
        end
