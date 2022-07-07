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
c       BS              Base subroutine to use (real*8 external function)
c       BC		Base coeff such that data= SUM{i} BC(i)*BS(i)
c       BP              Base parameters for spatial and temporal positions
c  Output:
c	ddat		ddat value = SUM{i} DW(i)
c       DW              WORK ARRAY that contains BS output
c                       (i.e. DW(i)=BC(i)*BS(i) most of the time)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	recursive subroutine cpt_dat(nub,nb,BS,BC,BP,ddat,DW)
c
        implicit none
c
        integer :: nub,nb,ib
        real*8 :: BC(*),BP(*),ddat,DW(*)
c
        external BS
c
        call BS('f',nub,nb,BC,BP,DW) 
c
        ddat=0.0d0
        do ib=1,nb
          ddat=ddat+DW(ib)
        enddo
c
        return
        end
