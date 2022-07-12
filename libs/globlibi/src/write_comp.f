ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine write_comp
c		V. Lesur	07/03/2002
c
c       write a file of nl data line.
c       The line before the first data line will be '##########'
c       previous lines must start with a '#' and may contain comments
c       data lines have nd+3 real records:
c       point position (nd real), data, fit, (data-fit) 
c       The data point printed are selected by type
c
c       input:
c         nom		file name to write
c         nd            number of space dimmension
c         nmin          starting count
c         nmax          finishing count
c         ddat(nd+1,*)  position and data (dim (nd+1,*))      
c         xyzf(*)       fit to data (dim nl)
c         diff(*)       two sums over data points
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine write_comp(nom,nd,nmin,nmax,ddat,xyzf,diff)
c
        implicit none
c
        include 'mpif.h'
c
	character nom*100,theformat*15
        integer i,j,nmin,nmax,nd,nd1,nd2
        real*8 ddat(*),xyzf(*),diff(*)
        integer rank, ierr
c
        call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
c
        if (rank .eq. 0) then
          open(10,file=nom,status='unknown')
          write(10,'(A)')'##########'
        else
	  open(10,file=nom,position='append',status='unknown')
        endif
c
        diff(1)=0.0d0
        diff(2)=0.0d0
        nd1=int((nd+3)/10.)
        nd2=(nd+3)-nd1*10
        theformat='('//char(nd2+48)//'e17.9)'
        if(nd1.ne.0)theformat='('//char(nd1+48)//char(nd2+48)//'e17.9)'
        nd1=nd+1
c
        do i=nmin,nmax
          write(10,theformat)(ddat(j+(i-1)*nd1),j=1,nd1)
     >              ,xyzf(i),ddat(i*nd1)-xyzf(i)
          diff(1)=diff(1)+ddat(i*nd1)-xyzf(i)
          diff(2)=diff(2)+(ddat(i*nd1)-xyzf(i))**2
        enddo
c
        close(10)
c
        return
        end