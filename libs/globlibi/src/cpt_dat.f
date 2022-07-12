cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       subroutine cpt_dat
c               V. Lesur 09/02/2005
c
c   18.11.2014: DATA variable changed to ddat to avoid possible 
c               confusion by compilers
c
c       recusive subroutine, should be compiled using F90
c       previously fmxyz.f.
c
c       Forward modelling for a single data value
c       Output also one row of the condition matrix
c
c  Input:
c       nub             base number 
c       nb              Number of bases
c       _base           Base subroutine to use (real*8 external function)
c       bc		Base coeff such that data= SUM{i} BC(i)*BS(i)
c       bp              Base parameters for spatial and temporal positions
c  Output:
c	ddat		ddat value = SUM{i} DW(i)
c       dw              WORK ARRAY that contains BS output
c                       (i.e. DW(i)=BC(i)*BS(i) most of the time)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine cpt_dat(nub,nb,sub_base,bc,bp,ddat,dw)
c
        implicit none
c
        integer :: nub,nb,ib
        real*8 :: bc(*),bp(*),ddat,dw(*)
c
        external sub_base
c
        call sub_base('f',nub,nb,bc,bp,dw) 
c
        ddat=0.0d0
        do ib=1,nb
          ddat=ddat+dw(ib)
        enddo
c
        return
        end